#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"
#include "../src/kmer.h"
#include "../src/pattern.h"
#include "../src/transbuilder.h"

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

  try {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "display this help message")
      ("kmerlen,k", po::value<int>()->default_value(15), "length of k-mers in de Bruijn graph")
      ("tandem,t", po::value<int>(), "reject local tandem duplications & inverted repeats up to this length")
      ("invrep,i", po::value<int>()->default_value(5), "reject nonlocal inverted repeats of this length (separated by at least 2 bases)")
      ("exclude,x", po::value<vector<string> >(), "motif(s) to exclude")
      ("source,s", po::value<vector<string> >(), "source motif(s): machine can start in this state, but will never enter it")
      ("control,c", po::value<int>()->default_value(0), "number of control words")
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
    
    const Pos len = vm.at("kmerlen").as<int>();
    Assert (len <= 31, "Maximum context is 31 bases");

    TransBuilder builder (len);

    if (vm.count("tandem"))
      builder.maxTandemRepeatLen = vm.at("tandem").as<int>();
    builder.invertedRepeatLen = vm.at("invrep").as<int>();

    getMotifs (vm, "exclude", builder.excludedMotif, builder.excludedMotifRevComp);
    getMotifs (vm, "source", builder.sourceMotif, builder.excludedMotifRevComp);

    builder.controlWords = vm.at("control").as<int>();
    
    // build transducer
    builder.removeRepeats();
    builder.pruneDeadEnds();
    builder.pruneUnreachable();
    builder.indexStates();

    // Output the transducer
    builder.output (cout);
    
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}
