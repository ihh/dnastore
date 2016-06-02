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
#include "../src/mutator.h"
#include "../src/fwdback.h"
#include "../src/viterbi.h"

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
      ("elim-trans", "eliminate degenerate transitions")
      ("controls,c", po::value<int>()->default_value(4), "number of control words")
      ("print-controls", "print control words")
      ("no-start", "do not use a control word at start of encoded sequence")
      ("no-end", "do not use a control word at end of encoded sequence")
      ("delay,y", "build delayed machine")
      ("rate,R", "calculate compression rate")
      ("dot", "print in Graphviz format")
      ("token-info", "print descriptions of input tokens")
      ("load-machine,L", po::value<string>(), "load machine from JSON file")
      ("save-machine,S", po::value<string>(), "save machine to JSON file")
      ("compose-machine,C", po::value<vector<string> >(), "load machine from JSON file and compose in front of primary machine")
      ("encode-file,e", po::value<string>(), "encode binary file to FASTA on stdout")
      ("decode-file,d", po::value<string>(), "decode FASTA file to binary on stdout")
      ("encode-string,E", po::value<string>(), "encode ASCII string to FASTA on stdout")
      ("decode-string,D", po::value<string>(), "decode DNA sequence to binary on stdout")
      ("encode-bits,b", po::value<string>(), "encode string of bits and control symbols to FASTA on stdout")
      ("decode-bits,B", po::value<string>(), "decode DNA sequence to string of bits and control symbols on stdout")
      ("decode-viterbi,V", po::value<string>(), "decode FASTA file using Viterbi algorithm")
      ("raw,r", "strip headers from FASTA output; just print raw sequence")
      ("error-sub-prob", po::value<double>()->default_value(.1), "substitution probability for error model")
      ("error-iv-ratio", po::value<double>()->default_value(10), "transition/transversion ratio for error model")
      ("error-dup-prob", po::value<double>()->default_value(.01), "tandem duplication probability for error model")
      ("error-del-open", po::value<double>()->default_value(.01), "deletion opening probability for error model")
      ("error-del-ext", po::value<double>()->default_value(.01), "deletion extension probability for error model")
      ("error-global", "force global alignment in error model (disallow partial reads)")
      ("error-file,F", po::value<string>(), "load error model from file")
      ("fit-error,f", po::value<string>(), "train error model on Stockholm database of pairwise alignments and print to stdout")
      ("error-counts", po::value<string>(), "estimate posterior expected counts of various different types of error from Stockholm database")
      ("strict-guides", "treat alignments in Stockholm database as strict truth, not just hints")
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

    if (vm.count("elim-trans"))
      builder.keepDegenerates = false;

    builder.nControlWords = vm.at("controls").as<int>();
    builder.controlWordAtStart = !vm.count("no-start");
    builder.controlWordAtEnd = !vm.count("no-end");
    builder.buildDelayedMachine = vm.count("delay");

    MutatorParams mut;
    if (vm.count("error-file")) {
      mut = MutatorParams::fromFile (vm.at("error-file").as<string>().c_str());
      LogThisAt(6,"Loaded error model:\n" << mut.asJSON());
    } else {
      mut.initMaxDupLen (len / 2);
      mut.pTanDup = vm.at("error-dup-prob").as<double>();
      mut.pDelOpen = vm.at("error-del-open").as<double>();
      mut.pDelExtend = vm.at("error-del-ext").as<double>();
      const double subProb = vm.at("error-sub-prob").as<double>();
      const double ivRatio = vm.at("error-iv-ratio").as<double>();
      mut.pTransition = subProb * ivRatio / (1 + ivRatio);
      mut.pTransversion = subProb / (1 + ivRatio);
      mut.local = vm.count("error-global") ? false : true;
      LogThisAt(6,"Command line-specified error model:\n" << mut.asJSON());
    }

    const bool rawSeqOutput = vm.count("raw");
    const bool strictAlignments = vm.count("strict-guides");
    
    if (vm.count("fit-error")) {
      const list<Stockholm> db = readStockholmDatabase (vm.at("fit-error").as<string>().c_str());
      MutatorCounts prior (mut);
      prior.initLaplace();
      const MutatorParams fitMut = baumWelchParams (mut, prior, db, strictAlignments);
      fitMut.writeJSON (cout);

    } else if (vm.count("error-counts")) {
      const list<Stockholm> db = readStockholmDatabase (vm.at("error-counts").as<string>().c_str());
      LogProb ll;
      const MutatorCounts counts = expectedCounts (mut, db, ll, strictAlignments);
      counts.writeJSON (cout);

    } else {
      // build, or load, transducer
      const bool loadMachine = vm.count("load-machine");
      Machine machine = loadMachine
	? Machine::fromFile(vm.at("load-machine").as<string>().c_str())
	: builder.makeMachine();

      if (!loadMachine && vm.count("print-controls"))
	cout << "Control words: " << join(builder.controlWordString) << endl;

      // pre-compose transducers
      if (vm.count("compose-machine")) {
	const vector<string> comps = vm.at("compose-machine").as<vector<string> >();
	for (auto iter = comps.rbegin(); iter != comps.rend(); ++iter) {
	  LogThisAt(3,"Pre-composing with " << *iter << endl);
	  machine = Machine::compose (Machine::fromFile((*iter).c_str()), machine);
	}
      }

      // save transducer
      if (vm.count("save-machine")) {
	const string savefile = vm.at("save-machine").as<string>();
	if (savefile == "-")
	  machine.writeJSON (cout);
	else {
	  ofstream out (savefile);
	  machine.writeJSON (out);
	}
      }

      // encoding or decoding?
      if (vm.count("encode-file")) {
	const string filename = vm.at("encode-file").as<string>();
	ifstream infile (filename, std::ios::binary);
	if (!infile)
	  throw runtime_error ("Binary file not found");
	FastaWriter writer (cout, rawSeqOutput ? NULL : filename.c_str());
	Encoder<FastaWriter> encoder (machine, writer);
	encoder.encodeStream (infile);
	
      } else if (vm.count("decode-file")) {
	const vguard<FastSeq> fastSeqs = readFastSeqs (vm.at("decode-file").as<string>().c_str());
	BinaryWriter writer (cout);
	Decoder<BinaryWriter> decoder (machine, writer);
	for (auto& fs: fastSeqs)
	  decoder.decodeString (fs.seq);

      } else if (vm.count("encode-string")) {
	FastaWriter writer (cout, rawSeqOutput ? NULL : "ASCII_string");
	Encoder<FastaWriter> encoder (machine, writer);
	encoder.encodeString (vm.at("encode-string").as<string>());
      
      } else if (vm.count("decode-string")) {
	BinaryWriter writer (cout);
	Decoder<BinaryWriter> decoder (machine, writer);
	decoder.decodeString (vm.at("decode-string").as<string>());

      } else if (vm.count("encode-bits")) {
	FastaWriter writer (cout, rawSeqOutput ? NULL : "bit_string");
	Encoder<FastaWriter> encoder (machine, writer);
	encoder.encodeSymbolString (vm.at("encode-bits").as<string>());
      
      } else if (vm.count("decode-bits")) {
	Decoder<ostream> decoder (machine, cout);
	decoder.decodeString (vm.at("decode-bits").as<string>());
	decoder.close();

	cout << endl;

      } else if (vm.count("decode-viterbi")) {
	const auto decoded = decodeFastSeqs (vm.at("decode-viterbi").as<string>().c_str(), machine, mut);
	if (rawSeqOutput)
	  for (const auto& fs: decoded)
	    cout << fs.seq << endl;
	else
	  writeFastaSeqs (cout, decoded);
	
      } else if (vm.count("rate")) {
	// Output statistics
	const auto charBases = machine.expectedBasesPerInputSymbol("01$");
	vguard<string> cbstr;
	for (const auto& cb: charBases)
	  cbstr.push_back (Machine::charToString(cb.first) + ": " + to_string(cb.second));
	cout << "Expected bases/symbol: { " << join(cbstr,", ") << " }" << endl;

      } else {
	// Output the transducer
	if (vm.count("dot"))
	  machine.writeDot (cout);
	else if (vm.count("token-info"))
	  cout << machine.inputDescriptionTable();
	else if (!vm.count("save-machine"))
	  machine.write (cout);
      }
    }

#ifndef DEBUG
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
#endif /* DEBUG */
  
  return EXIT_SUCCESS;
}
