#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"

using namespace std;

namespace po = boost::program_options;

typedef unsigned long long Kmer;
typedef unsigned char Base;
typedef int Pos;

typedef unsigned char EdgeFlags;
typedef unsigned long long State;

const char* alphabet = "AGTC";
inline char baseToChar (Base base) {
  return alphabet[base & 3];
}

inline Base charToBase (char c) {
  const char* s = strchr (alphabet, toupper(c));
  Require (s != NULL, "%c is not a nucleotide character", c);
  return (Base) (s - alphabet);
}

inline Base getBase (Kmer kmer, Pos pos) {
  return (kmer >> ((pos - 1) << 1)) & 3;
}

inline Kmer setBase (Kmer kmer, Pos pos, Base base) {
  const int shift = (pos - 1) << 1;
  return (kmer & (((Kmer) -1) ^ (3 << shift))) | (((Kmer) base) << shift);
}

inline Base complementBase (Base b) {
  return b ^ 2;
}

inline Kmer makeTransition (Kmer kmer, Pos pos) {
  return kmer ^ (1 << ((pos - 1) << 1));
}

inline string kmerString (Kmer kmer, Pos len) {
  string s (len, '*');
  for (Pos i = 1; i <= len; ++i)
    s[len-i] = baseToChar(getBase(kmer,i));
  return s;
}

inline Kmer stringToKmer (const string& s) {
  Kmer kmer = 0;
  for (Pos i = 0; i <= s.size(); ++i)
    kmer = setBase (kmer, s.size() - i, charToBase (s[i]));
  return kmer;
}

inline bool isTransition (Base x, Base y) {
  return x != y && (x & 2) == (y & 2);
}

inline bool isTransversion (Base x, Base y) {
  return x != y && (x & 2) != (y & 2);
}

inline bool isComplement (Base x, Base y) {
  return y == complementBase(x);
}

inline bool isGC (Base x) {
  return (x & 1) == 1;
}

inline double gcContent (Kmer kmer, Pos len) {
  int gc = 0;
  for (Pos pos = 1; pos <= len; ++pos)
    if (isGC (getBase (kmer, pos)))
      ++gc;
  return gc / (double) len;
}

inline double gcNonuniformity (Kmer kmer, Pos len) {
  const double gc = gcContent (kmer, len);
  return abs (gc - .5);
}

inline double kmerEntropy (Kmer kmer, Pos len) {
  vguard<int> freq (4);
  for (Pos pos = 1; pos <= len; ++pos)
    ++freq[getBase(kmer,pos)];
  double S = 0;
  for (Base b = 0; b < 3; ++b)
    if (freq[b])
      S -= freq[b] * log (freq[b] / (double) len);
  return S / log(2);
}

inline bool kmerEqualOrBetter (Kmer x, Kmer y, Pos len) {
  const double xgc = gcNonuniformity(x,len);
  const double ygc = gcNonuniformity(y,len);
  return xgc == ygc ? (kmerEntropy(x,len) >= kmerEntropy(y,len)) : (xgc < ygc);
}

inline Kmer kmerMask (Pos len) {
  // 4^len - 1 = 2^(2*len) - 1 = (1 << (2*len)) - 1 = (1 << (len << 1)) - 1
  return (((Kmer) 1) << (len << 1)) - 1;
}

inline Kmer kmerSub (Kmer kmer, Pos start, Pos len) {
  return (kmer >> ((start - 1) << 1)) & kmerMask (len);
}

inline string kmerSubstring (Kmer kmer, Pos start, Pos len) {
  return kmerString (kmerSub (kmer, start, len), len);
}

inline int kmerLeftCoord (Pos pos, Pos len) {
  return len - pos + 1;
}

inline string kmerSubCoords (Pos start, Pos len, Pos kmerLen) {
  return string("[") + to_string(kmerLeftCoord(start+len-1,kmerLen)) + (len > 1 ? (string("..") + to_string(kmerLeftCoord(start,kmerLen))) : string()) + "]";
}

inline string kmerSubAt (Kmer kmer, Pos start, Pos len, Pos kmerLen) {
  return kmerString (kmerSub (kmer, start, len), len) + kmerSubCoords(start,len,kmerLen);
}

inline Kmer kmerRevComp (Kmer kmer, Pos len) {
  Kmer rc = 0;
  for (Pos i = 1; i <= len; ++i)
    rc = (rc << 2) | complementBase (getBase(kmer,i));
  return rc;
}

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

inline EdgeFlags edgeFlags (Kmer kmer, Pos len, vguard<Kmer>& kmers, const vguard<bool>& validFlag) {
  EdgeFlags f = 0;
  for (size_t n = 0; n < 4; ++n)
    if (validFlag[kmers[n]])
      f = f | (1 << n);
  return f;
}

inline EdgeFlags outgoingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& outgoing, const vguard<bool>& validFlag) {
  getOutgoing (kmer, len, outgoing);
  return edgeFlags (kmer, len, outgoing, validFlag);
}

inline EdgeFlags incomingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& incoming, const vguard<bool>& validFlag) {
  getIncoming (kmer, len, incoming);
  return edgeFlags (kmer, len, incoming, validFlag);
}

int edgeFlagsToCountLookup[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
inline int edgeFlagsToCount (EdgeFlags flags) {
  return edgeFlagsToCountLookup[flags & 0xf];
}

inline void pruneDeadEnds (Kmer kmer, Pos len, vguard<Kmer>& in, vguard<Kmer>& out, vguard<bool>& valid) {
  if (valid[kmer]) {
    const EdgeFlags inFlags = incomingEdgeFlags(kmer,len,in,valid);
    const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,valid);
    if (inFlags == 0 || outFlags == 0) {
      LogThisAt(5,"Pruning " << kmerString(kmer,len) << endl);
      valid[kmer] = 0;
      for (auto kmerIn: in)
	pruneDeadEnds (kmerIn, len, in, out, valid);
      for (auto kmerOut: out)
	pruneDeadEnds (kmerOut, len, in, out, valid);
    } else
      LogThisAt(9,"Keeping " << kmerString(kmer,len) << " with " << edgeFlagsToCount(inFlags) << " incoming and " << edgeFlagsToCount(outFlags) << " outgoing edges" << endl);
  }
}

inline list<Kmer> pruneDeadEnds (const list<Kmer>& kmers, Pos len, vguard<bool>& valid) {
  const Kmer maxKmer = kmerMask(len);
  vguard<Kmer> incoming (4), outgoing (4);
  ProgressLog (plogPrune, 1);
  plogPrune.initProgress ("Pruning dead ends");
  for (auto kmer: kmers)
    pruneDeadEnds (kmer, len, incoming, outgoing, valid);
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

inline bool containsMotif (Kmer seq, Pos len, Kmer motif, Pos motifLen, const char* desc) {
  for (Pos i = len - motifLen + 1; i >= 1; --i)
    if (kmerSub (seq, i, motifLen) == motif) {
      LogThisAt(4,"Rejecting " << kmerString(seq,len) << " because it contains " << kmerSubAt(seq,i,motifLen,len) << " (" << desc << ")" << endl);
      return true;
    }

  return false;
}

inline bool containsMotifs (Kmer seq, Pos len, const vguard<Kmer>& motifs, const vguard<Pos>& motifLengths, const char* desc) {
  for (size_t n = 0; n < motifs.size(); ++n)
    if (containsMotif (seq, len, motifs[n], motifLengths[n], desc))
      return true;
  return false;
}

inline bool hasExactTandemRepeat (Kmer seq, Pos len, Pos maxRepeatLen) {
  for (Pos repeatLen = 1; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = len - 2*repeatLen + 1; i >= 1; --i)
      if (kmerSub (seq, i, repeatLen) == kmerSub (seq, i + repeatLen, repeatLen)) {
	const int logLevel = max(5,8-repeatLen);
	LogThisAt(logLevel,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,i+repeatLen,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (" << (repeatLen == 1 ? "repeated base" : "exact tandem repeat") << ")" << endl);
	return true;
      }
  return false;
}

inline bool hasExactLocalInvertedRepeat (Kmer seq, Pos len, Pos minRepeatLen, Pos maxRepeatLen) {
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos repeatLen = minRepeatLen; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = len - 2*repeatLen + 1; i >= 1; --i) {
      const Kmer invRep = kmerSub (rc, len - i + 1, repeatLen);
      if (invRep == kmerSub (seq, i + repeatLen, repeatLen)) {
	const int logLevel = max(5,8-repeatLen);
	LogThisAt(logLevel,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,i+repeatLen,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (palindrome)" << endl);
	return true;
      }
    }
  return false;
}

inline bool hasExactNonlocalInvertedRepeat (Kmer seq, Pos len, Pos repeatLen, Pos minSeparation) {
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos i = len - repeatLen*2 - minSeparation; i > 0; --i) {
    const Kmer invRep = kmerSub (rc, len - i + 1, repeatLen);
    const Pos jMin = i + repeatLen + minSeparation;
    for (Pos j = len - repeatLen + 1; j >= jMin; --j)
      if (invRep == kmerSub (seq, j, repeatLen)) {
	LogThisAt(4,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,j,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (exact inverted repeat)" << endl);
	return true;
      }
  }
  return false;
}

int main (int argc, char** argv) {

  try {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "display this help message")
      ("kmerlen,k", po::value<int>()->default_value(15), "length of k-mers in de Bruijn graph")
      ("tandem,t", po::value<int>()->default_value(5), "reject local tandem duplications & inverted repeats up to this length")
      ("invrep,i", po::value<int>()->default_value(5), "reject nonlocal inverted repeats of this length (separated by at least 2 bases)")
      ("exclude,x", po::value<vector<string> >(), "motif(s) to exclude")
      ("verbose,v", po::value<int>()->default_value(2), "verbosity level")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    logger.setVerbose (vm.at("verbose").as<int>());
    
    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }

    const Pos len = vm.at("kmerlen").as<int>();
    Assert (len <= 31, "Maximum context is 31 bases");
    
    const Pos maxTandemRepeatLen = vm.at("tandem").as<int>();
    const Pos invertedRepeatLen = vm.at("invrep").as<int>();

    vguard<Kmer> excludedMotif, excludedMotifRevComp;
    vguard<Pos> excludedMotifLen;
    if (vm.count("exclude"))
      for (const auto& x: vm.at("exclude").as<vector<string> >()) {
	const Kmer motif = stringToKmer (x);
	const Pos motifLen = x.size();
	excludedMotif.push_back (motif);
	excludedMotifRevComp.push_back (kmerRevComp (motif, motifLen));
	excludedMotifLen.push_back (motifLen);
      }

    const Kmer maxKmer = kmerMask(len);

    // Remove repeats
    ProgressLog (plogReps, 1);
    plogReps.initProgress ("Filtering %d-mer repeats", len);
    vguard<bool> kmerValid (maxKmer+1);
    unsigned long long nKmersWithoutReps = 0;
    list<Kmer> kmersWithoutReps;
    for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
      plogReps.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      
      if (!containsMotifs(kmer,len,excludedMotif,excludedMotifLen,"excluded motif")
	  && !containsMotifs(kmer,len,excludedMotifRevComp,excludedMotifLen,"revcomp of excluded motif")
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
    const list<Kmer> prunedKmers = pruneDeadEnds (kmersWithoutTransitions, len, kmerValid);
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
      const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,kmerValid);
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
