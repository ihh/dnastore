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
      ("source,s", po::value<vector<string> >(), "source motif(s): machine can start in this state, but will never enter it")
      ("keepdeg,k", "keep degenerate transitions")
      ("control,c", po::value<int>()->default_value(4), "number of control words")
      ("encode-file,e", po::value<string>(), "encode binary file to FASTA on stdout")
      ("decode-file,d", po::value<string>(), "decode FASTA file to binary on stdout")
      ("encode-string,E", po::value<string>(), "encode ASCII string to FASTA on stdout")
      ("decode-string,D", po::value<string>(), "decode DNA sequence to binary on stdout")
      ("encode-bits,b", po::value<string>(), "encode bit-string to FASTA on stdout")
      ("decode-bits,B", po::value<string>(), "decode DNA sequence to bit-string on stdout")
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

    if (vm.count("keepdeg"))
      builder.keepDegenerates = true;
    builder.nControlWords = vm.at("control").as<int>();
    
    // build transducer
    const Machine machine = builder.makeMachine();
    
    // encoding or decoding?
    if (vm.count("encode-file")) {
      ifstream infile (vm.at("encode-file").as<string>(), std::ios::binary);
      if (!infile)
	throw runtime_error ("Binary file not found");
      Encoder encoder (machine);
      FastaWriter writer (cout);
      encoder.encodeStream (infile, writer);
      
    } else if (vm.count("decode-file")) {
      const vguard<FastSeq> fastSeqs = readFastSeqs (vm.at("decode").as<string>().c_str());
      Decoder decoder (machine);
      BinaryWriter writer (cout);
      for (auto& fs: fastSeqs)
	decoder.decodeString (fs.seq, writer);

    } else if (vm.count("encode-string")) {
      Encoder encoder (machine);
      FastaWriter writer (cout);
      encoder.encodeString (vm.at("encode-string").as<string>(), writer);
      
    } else if (vm.count("decode-string")) {
      Decoder decoder (machine);
      BinaryWriter writer (cout);
      decoder.decodeString (vm.at("decode-string").as<string>(), writer);

    } else if (vm.count("encode-bits")) {
      Encoder encoder (machine);
      FastaWriter writer (cout);
      encoder.encodeSymbolString (vm.at("encode-bits").as<string>(), writer);
      
    } else if (vm.count("decode-bits")) {
      Decoder decoder (machine);
      decoder.decodeString (vm.at("decode-bits").as<string>(), cout);
      cout << endl;

    } else {
      // Output the transducer
      machine.write (cout);

      // Output statistics
      const double basesPerBit = machine.expectedBasesPerBit();
      LogThisAt(1,"Expected bases/bit: " << basesPerBit << endl);
      if (builder.nControlWords)
	LogThisAt(1,"Expected bases/control: " << builder.expectedBasesPerControlChar() << endl);
    }

#ifndef DEBUG
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
#endif /* DEBUG */
  
  return EXIT_SUCCESS;
}
