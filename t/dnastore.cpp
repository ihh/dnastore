#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"
#include "../src/kmer.h"
#include "../src/pattern.h"

using namespace std;

namespace po = boost::program_options;

inline void getOutgoing (Kmer kmer, Pos len, vguard<Kmer>& outgoing) {
  Assert (outgoing.size() == 4, "oops");
  const Kmer prefix = (kmer << 2) & kmerMask(len);
  iota (outgoing.begin(), outgoing.end(), prefix);
}

inline void getIncoming (Kmer kmer, Pos len, vguard<Kmer>& incoming) {
  Assert (incoming.size() == 4, "oops");
  const Kmer prefix = kmer >> 2;
  const int shift = (len - 1) << 1;
  for (Base b = 0; b < 4; ++b)
    incoming[b] = prefix | (((Kmer) b) << shift);
}

inline EdgeFlags outgoingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& outgoing, const vguard<bool>& validFlag, const set<KmerLen>& sources) {
  getOutgoing (kmer, len, outgoing);
  EdgeFlags f = 0;
  for (size_t n = 0; n < 4; ++n)
    if (validFlag[outgoing[n]] && !endsWithMotif(outgoing[n],len,sources))
      f = f | (1 << n);
  return f;
}

inline EdgeFlags incomingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& incoming, const vguard<bool>& validFlag) {
  getIncoming (kmer, len, incoming);
  EdgeFlags f = 0;
  for (size_t n = 0; n < 4; ++n)
    if (validFlag[incoming[n]])
      f = f | (1 << n);
  return f;
}

int edgeFlagsToCountLookup[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
inline int edgeFlagsToCount (EdgeFlags flags) {
  return edgeFlagsToCountLookup[flags & 0xf];
}

inline void pruneDeadEnds (Kmer kmer, Pos len, vguard<bool>& valid, const set<KmerLen>& sources) {
  vguard<Kmer> in (4), out (4);
  if (valid[kmer] && !endsWithMotif(kmer,len,sources)) {
    const EdgeFlags inFlags = incomingEdgeFlags(kmer,len,in,valid);
    const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,valid,sources);
    const bool prune = inFlags == 0 || outFlags == 0;
    LogThisAt(9,(prune ? "Pruning" : "Keeping") << " " << kmerString(kmer,len) << " with " << edgeFlagsToCount(inFlags) << " incoming and " << edgeFlagsToCount(outFlags) << " outgoing edges" << endl);
    if (prune) {
      valid[kmer] = 0;
      for (auto kmerIn: in)
	pruneDeadEnds (kmerIn, len, valid, sources);
      for (auto kmerOut: out)
	pruneDeadEnds (kmerOut, len, valid, sources);
    }
  }
}

inline list<Kmer> pruneDeadEnds (const list<Kmer>& kmers, Pos len, vguard<bool>& valid, const set<KmerLen>& sources) {
  const Kmer maxKmer = kmerMask(len);
  ProgressLog (plogPrune, 1);
  plogPrune.initProgress ("Pruning dead ends");
  for (auto kmer: kmers)
    pruneDeadEnds (kmer, len, valid, sources);
  const unsigned long long nKmers = kmers.size();
  unsigned long long nPruned = 0, nUnpruned = 0;
  list<Kmer> unprunedKmers;
  for (auto kmer: kmers) {
    ++nPruned;
    plogPrune.logProgress (nPruned / (double) nKmers, "sequence %llu/%llu", nPruned, nKmers);
    if (valid[kmer]) {
      unprunedKmers.push_back (kmer);
      ++nUnpruned;
    }
  }
  LogThisAt(2,"Dead-end pruning left " << nUnpruned << " " << len << "-mers (" << setprecision(2) << 100*(double)nUnpruned/(1.+(double)maxKmer) << "%)" << endl);
  return unprunedKmers;
}

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
      ("keep,k", po::value<bool>(), "keep transition degeneracies")
      ("verbose,v", po::value<int>()->default_value(2), "verbosity level")
      ("log", po::value<vector<string> >(), "log everything in this function")
      ("nocolor", po::value<bool>(), "log in monochrome")
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
    
    const Pos maxTandemRepeatLen = vm.count("tandem") ? vm.at("tandem").as<int>() : (len / 2);
    const Pos invertedRepeatLen = vm.at("invrep").as<int>();

    set<KmerLen> excludedMotif, excludedMotifRevComp;
    getMotifs (vm, "exclude", excludedMotif, excludedMotifRevComp);

    set<KmerLen> sourceMotif;
    getMotifs (vm, "source", sourceMotif, excludedMotifRevComp);

    const bool keepDegeneracies = vm.count("keep");
    
    // Remove repeats
    const Kmer maxKmer = kmerMask(len);
    vguard<bool> kmerValid (maxKmer+1);
    unsigned long long nKmersWithoutReps = 0;
    list<Kmer> kmersWithoutReps;
    ProgressLog (plogReps, 1);
    plogReps.initProgress ("Filtering %d-mer repeats", len);
    for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
      plogReps.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      
      if (!endsWithMotif(kmer,len,excludedMotif,"excluded motif")
	  && !endsWithMotif(kmer,len,excludedMotifRevComp,"revcomp of excluded motif")
	  && !hasExactTandemRepeat(kmer,len,maxTandemRepeatLen)
	  && !hasExactLocalInvertedRepeat(kmer,len,3,maxTandemRepeatLen)
	  && !hasExactNonlocalInvertedRepeat(kmer,len,invertedRepeatLen,2)) {
	LogThisAt(9,"Accepting " << kmerString(kmer,len) << endl);
	kmerValid[kmer] = true;
	kmersWithoutReps.push_back (kmer);
	++nKmersWithoutReps;
      }
    }
    LogThisAt(2,"Found " << nKmersWithoutReps << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutReps/(1.+(double)maxKmer) << "%)" << endl);

    // Eliminate transition redundancies where possible
    list<Kmer> kmersWithoutTransitions;
    for (Kmer kmer: kmersWithoutReps) {
      const Kmer transKmer = makeTransition(kmer,1);
      if (kmerValid[kmer] && kmerValid[transKmer] && kmerEqualOrBetter(transKmer,kmer,len)) {
	kmerValid[kmer] = false;
	LogThisAt(4,"Eliminating "
		  << kmerString(kmer,len) << " (" << setprecision(2) << 100*gcContent(kmer,len) << "% GC, " << setw(3) << kmerEntropy(kmer,len) << " bits)"
		  << " since we also have "
		  << kmerString(transKmer,len) << " (" << setprecision(2) << 100*gcContent(transKmer,len) << "% GC, " << setw(3) << kmerEntropy(transKmer,len) << " bits)"
		  << endl);
      } else
	kmersWithoutTransitions.push_back (kmer);
    }
    const unsigned long long nKmersWithoutTransitions = kmersWithoutTransitions.size();
    LogThisAt(2,"Removed " << (nKmersWithoutReps - nKmersWithoutTransitions) << " transition redundancies leaving " << nKmersWithoutTransitions << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutTransitions/(1.+(double)maxKmer) << "%)" << endl);
    
    // Remove dead ends
    const list<Kmer> prunedKmers = pruneDeadEnds (kmersWithoutTransitions, len, kmerValid, sourceMotif);
    const unsigned long long nPrunedKmers = prunedKmers.size();

    Require (nPrunedKmers > 0, "No valid %d-mers left!", len);

    // Index the states
    map<Kmer,State> kmerState;
    State n = 0;
    for (auto kmer: prunedKmers)
      kmerState[kmer] = ++n;

    // Output the transducer
    vguard<Kmer> out (4);
    vguard<char> outChar;
    vguard<State> outState;
    for (auto kmer: prunedKmers) {
      const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,kmerValid,sourceMotif);
      outChar.clear();
      outState.clear();
      for (size_t n = 0; n < 4; ++n)
	if (outFlags & (1 << n)) {
	  outChar.push_back (baseToChar(n));
	  outState.push_back (kmerState.at(out[n]));
	}
      cout << "#" << kmerState.at(kmer) << " " << kmerString(kmer,len);
      if (outChar.size() == 1)
	cout << " /" << outChar[0] << "->#" << outState[0];
      else if (outChar.size() == 2)
	cout << " 0/" << outChar[0] << "->#" << outState[0]
	     << " 1/" << outChar[1] << "->#" << outState[1];
      else if (outChar.size() == 3)
	cout << " 00/" << outChar[0] << "->#" << outState[0]
	     << " 01/" << outChar[1] << "->#" << outState[1]
	     << " 1/" << outChar[2] << "->#" << outState[2];
      cout << endl;
    }

  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}
