#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"
#include "../src/kmer.h"
#include "../src/pattern.h"
#include "../src/builder.h"
#include "../src/encoder.h"
#include "../src/decoder.h"
#include "../src/fastseq.h"

using namespace std;

namespace po = boost::program_options;

void getMotifs (po::variables_map& vm, const char* arg, set<KmerLen>& motifs, set<KmerLen>& motifRevComps) {
  if (vm.count(arg))
    for (const auto& x: vm.at(arg).as<vector<string> >()) {
      const Kmer motif = stringToKmer (x);
      const Pos motifLen = x.size();
      motifs.insert (KmerLen (motif, motifLen));
      motifRevComps.insert (KmerLen (kmerRevComp (motif, motifLen), motifLen));
    }
}

int main (int argc, char** argv) {

#ifndef DEBUG
  try {
#endif /* DEBUG */
    
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "display this help message")
      ("length,l", po::value<int>()->default_value(12), "length of k-mers in de Bruijn graph")
      ("tandem,t", po::value<int>(), "reject local tandem duplications & inverted repeats up to this length")
      ("invrep,i", po::value<int>()->default_value(4), "reject nonlocal inverted repeats of this length (separated by at least 2 bases)")
      ("exclude,x", po::value<vector<string> >(), "motif(s) to exclude")
      ("source,o", po::value<vector<string> >(), "source motif(s): machine can start in this state, but will never enter it")
      ("keep-degen,k", "keep degenerate transitions")
      ("controls,c", po::value<int>()->default_value(4), "number of control words")
      ("no-start,s", "do not use a control word at start of encoded sequence")
      ("no-eof,f", "do not use a control word at end of encoded sequence")
      ("delay,y", "build delayed machine")
      ("load-machine,L", po::value<string>(), "load machine from JSON file")
      ("save-machine,S", po::value<string>(), "save machine to JSON file")
      ("encode-file,e", po::value<string>(), "encode binary file to FASTA on stdout")
      ("decode-file,d", po::value<string>(), "decode FASTA file to binary on stdout")
      ("encode-string,E", po::value<string>(), "encode ASCII string to FASTA on stdout")
      ("decode-string,D", po::value<string>(), "decode DNA sequence to binary on stdout")
      ("encode-bits,b", po::value<string>(), "encode string of bits and control symbols to FASTA on stdout")
      ("decode-bits,B", po::value<string>(), "decode DNA sequence to string of bits and control symbols on stdout")
      ("raw,r", "strip headers from FASTA output; just print raw sequence")
      ("verbose,v", po::value<int>()->default_value(2), "verbosity level")
      ("log", po::value<vector<string> >(), "log everything in this function")
      ("nocolor", "log in monochrome")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    // parse args
    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }

    logger.parseLogArgs (vm);
    
    const Pos len = vm.at("length").as<int>();
    Assert (len <= 31, "Maximum context is 31 bases");

    TransBuilder builder (len);

    if (vm.count("tandem"))
      builder.maxTandemRepeatLen = vm.at("tandem").as<int>();
    builder.invertedRepeatLen = vm.at("invrep").as<int>();

    getMotifs (vm, "exclude", builder.excludedMotif, builder.excludedMotifRevComp);
    getMotifs (vm, "source", builder.sourceMotif, builder.excludedMotifRevComp);

    if (vm.count("keep-degen"))
      builder.keepDegenerates = true;

    builder.nControlWords = vm.at("controls").as<int>();
    builder.controlWordAtStart = !vm.count("no-start");
    builder.controlWordAtEnd = !vm.count("no-eof");
    builder.buildDelayedMachine = vm.count("delay");

    // build, or load, transducer
    const Machine machine = vm.count("load-machine")
      ? Machine::fromFile(vm.at("load-machine").as<string>().c_str())
      : builder.makeMachine();
    
    // save transducer
    if (vm.count("save-machine")) {
      ofstream out (vm.at("save-machine").as<string>());
      machine.writeJSON (out);
    }
    
    // encoding or decoding?
    if (vm.count("encode-file")) {
      const string filename = vm.at("encode-file").as<string>();
      ifstream infile (filename, std::ios::binary);
      if (!infile)
	throw runtime_error ("Binary file not found");
      FastaWriter writer (cout, vm.count("raw") ? NULL : filename.c_str());
      Encoder<FastaWriter> encoder (machine, writer);
      encoder.encodeStream (infile);
	
    } else if (vm.count("decode-file")) {
      const vguard<FastSeq> fastSeqs = readFastSeqs (vm.at("decode-file").as<string>().c_str());
      BinaryWriter writer (cout);
      Decoder<BinaryWriter> decoder (machine, writer);
      for (auto& fs: fastSeqs)
	decoder.decodeString (fs.seq);

    } else if (vm.count("encode-string")) {
      FastaWriter writer (cout, vm.count("raw") ? NULL : "ASCII_string");
      Encoder<FastaWriter> encoder (machine, writer);
      encoder.encodeString (vm.at("encode-string").as<string>());
      
    } else if (vm.count("decode-string")) {
      BinaryWriter writer (cout);
      Decoder<BinaryWriter> decoder (machine, writer);
      decoder.decodeString (vm.at("decode-string").as<string>());

    } else if (vm.count("encode-bits")) {
      FastaWriter writer (cout, vm.count("raw") ? NULL : "bit_string");
      Encoder<FastaWriter> encoder (machine, writer);
      encoder.encodeSymbolString (vm.at("encode-bits").as<string>());
      
    } else if (vm.count("decode-bits")) {
      Decoder<ostream> decoder (machine, cout);
      decoder.decodeString (vm.at("decode-bits").as<string>());
      decoder.close();

      cout << endl;

    } else {
      // Output the transducer
      machine.write (cout);

      // Output statistics
      const auto charBases = machine.expectedBasesPerInputSymbol (true);
      vguard<string> cbstr;
      for (const auto& cb: charBases)
	cbstr.push_back (Machine::charToString(cb.first) + ": " + to_string(cb.second));
      LogThisAt(1,"Expected bases/symbol: { " << join(cbstr,", ") << " }" << endl);
    }

#ifndef DEBUG
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
#endif /* DEBUG */
  
  return EXIT_SUCCESS;
}
