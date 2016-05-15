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

const char* alphabet = "ATGC";
inline char baseToChar (Base base) {
  return alphabet[base & 3];
}

inline Base getBase (Kmer kmer, Pos pos) {
  return (kmer >> ((pos - 1) << 1)) & 3;
}

inline Kmer setBase (Kmer kmer, Pos pos, Base base) {
  const int shift = (pos - 1) << 1;
  return (kmer & (((Kmer) -1) ^ (3 << shift))) | (((Kmer) base) << shift);
}

inline string kmerString (Kmer kmer, Pos len) {
  string s (len, '*');
  for (Pos i = 1; i <= len; ++i)
    s[len-i] = baseToChar(getBase(kmer,i));
  return s;
}

inline bool isTransition (Base x, Base y) {
  return (x & 1) == (y & 1);
}

inline bool isComplement (Base x, Base y) {
  return (x ^ 3) == y;
}

inline bool isGC (Base x) {
  return (x & 2) == 1;
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
  return string("[") + to_string(kmerLeftCoord(start+len-1,kmerLen)) + ".." + to_string(kmerLeftCoord(start,kmerLen)) + "]";
}

inline string kmerSubAt (Kmer kmer, Pos start, Pos len, Pos kmerLen) {
  return kmerString (kmerSub (kmer, start, len), len) + kmerSubCoords(start,len,kmerLen);
}

inline Kmer kmerRevComp (Kmer kmer, Pos len) {
  Kmer rc = 0;
  for (Pos i = 1; i <= len; ++i)
    rc = (rc << 2) | getBase(kmer,i);
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

inline EdgeFlags outgoingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& outgoing, const vguard<bool>& validFlag) {
  getOutgoing (kmer, len, outgoing);
  EdgeFlags f = 0;
  for (size_t n = 0; n < 4; ++n)
    if (validFlag[outgoing[n]])
      f = f | (1 << n);
  return f;
}

int edgeFlagsToCountLookup[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
inline int edgeFlagsToCount (EdgeFlags flags) {
  return edgeFlagsToCountLookup[flags & 0xf];
}

void pruneKmer (Kmer kmer, Pos len, vguard<Kmer>& out, vguard<bool>& valid);
inline void pruneDeadEnds (Kmer kmer, Pos len, vguard<Kmer>& out, vguard<bool>& valid) {
  if (valid[kmer]) {
    const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,valid);
    if (outFlags == 0) {
      LogThisAt(9,"Pruning " << kmerString(kmer,len) << endl);
      pruneKmer(kmer,len,out,valid);
    } else
      LogThisAt(9,"Keeping " << kmerString(kmer,len) << " with " << edgeFlagsToCount(outFlags) << " outgoing edges" << endl);
  }
}

inline void pruneKmer (Kmer kmer, Pos len, vguard<Kmer>& out, vguard<bool>& valid) {
  valid[kmer] = 0;
  vguard<Kmer> incoming (4);
  getIncoming (kmer, len, incoming);
  for (auto kmerIn: incoming)
    pruneDeadEnds (kmerIn, len, out, valid);
}

inline list<Kmer> pruneDeadEnds (const list<Kmer>& kmers, Pos len, vguard<bool>& valid) {
  const Kmer maxKmer = kmerMask(len);
  vguard<Kmer> outgoing (4);
  ProgressLog (plogPrune, 1);
  plogPrune.initProgress ("Pruning dead ends");
  for (auto kmer: kmers)
    pruneDeadEnds (kmer, len, outgoing, valid);
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
  LogThisAt(2,"Dead-end pruning left " << nUnpruned << " candidate " << len << "-mers (" << setprecision(2) << 100*(double)nUnpruned/(1.+(double)maxKmer) << "%)" << endl);
  return unprunedKmers;
}

inline bool hasExactTandemRepeat (Kmer seq, Pos len, Pos maxRepeatLen) {
  for (Pos repeatLen = 1; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = len - repeatLen; i >= 1; --i)
      if (kmerSub (seq, i, repeatLen) == kmerSub (seq, i + repeatLen, repeatLen)) {
	const int logLevel = max(5,8-repeatLen);
	LogThisAt(logLevel,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,i+repeatLen,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (exact tandem repeat)" << endl);
	return true;
      }
  return false;
}

inline bool tooManyMismatches (Kmer seq1, Kmer seq2, Pos len, unsigned int maxMismatches) {
  unsigned int mismatches = 0;
  for (Pos i = 1; i < len; ++i)
    if (getBase(seq1,i) != getBase(seq2,i))
      if (++mismatches > maxMismatches)
	return true;
  return false;
}

inline bool hasExactInvertedRepeat (Kmer seq, Pos len, Pos repeatLen) {
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos i = len - repeatLen*2 - 2; i > 0; --i) {
    const Kmer invRep = kmerSub (rc, len - i + 1, repeatLen);
    const Pos jMin = i + repeatLen + 2;
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
      ("tandem,t", po::value<int>()->default_value(5), "reject local tandem duplications up to this length")
      ("invrep,i", po::value<int>()->default_value(5), "reject inverted repeats of this length (separated by at least 2 bases)")
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

    const Kmer maxKmer = kmerMask(len);

    // Remove repeats
    ProgressLog (plogReps, 1);
    plogReps.initProgress ("Filtering %d-mer repeats", len);
    vguard<bool> kmerValid (maxKmer+1);
    unsigned long long nKmersWithoutReps = 0;
    list<Kmer> kmersWithoutReps;
    for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
      plogReps.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      
      if (!hasExactTandemRepeat(kmer,len,maxTandemRepeatLen)
	  && !hasExactInvertedRepeat(kmer,len,invertedRepeatLen)) {
	LogThisAt(9,"Accepting " << kmerString(kmer,len) << endl);
	kmerValid[kmer] = true;
	kmersWithoutReps.push_back (kmer);
	++nKmersWithoutReps;
      }
    }
    LogThisAt(2,"Found " << nKmersWithoutReps << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutReps/(1.+(double)maxKmer) << "%)" << endl);

    // Remove dead ends
    const list<Kmer> prunedKmersWithoutReps = pruneDeadEnds (kmersWithoutReps, len, kmerValid);
    const unsigned long long nPrunedKmersWithoutReps = prunedKmersWithoutReps.size();

    Require (nPrunedKmersWithoutReps > 0, "No valid %d-mers left!", len);

    // Index the states
    map<Kmer,State> kmerState;
    State n = 0;
    for (auto kmer: prunedKmersWithoutReps)
      kmerState[kmer] = ++n;

    // Output the transducer
    vguard<Kmer> out (4);
    vguard<char> outChar;
    vguard<State> outState;
    for (auto kmer: prunedKmersWithoutReps) {
      const EdgeFlags outFlags = outgoingEdgeFlags(kmer,len,out,kmerValid);
      outChar.clear();
      outState.clear();
      for (size_t n = 0; n < 4; ++n)
	if (outFlags & (1 << n)) {
	  outChar.push_back (baseToChar(n));
	  outState.push_back (kmerState.at(out[n]));
	}
      cout << kmerState.at(kmer) << " " << kmerString(kmer,len);
      if (outChar.size() == 1)
	cout << " /" << outChar[0] << ":" << outState[0];
      else if (outChar.size() == 2)
	cout << " 0/" << outChar[0] << ":" << outState[0]
	     << " 1/" << outChar[1] << ":" << outState[1];
      else if (outChar.size() == 3)
	cout << " 00/" << outChar[0] << ":" << outState[0]
	     << " 01/" << outChar[1] << ":" << outState[1]
	     << " 1/" << outChar[2] << ":" << outState[2];
      cout << endl;
    }

  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}
