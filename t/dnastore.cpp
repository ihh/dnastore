#include <cstdlib>
#include <stdexcept>
#include <iomanip>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"

using namespace std;

namespace po = boost::program_options;

typedef unsigned long long Kmer;
typedef unsigned char Base;
typedef int Pos;

typedef unsigned char KmerFlags;
#define KmerValid 0x80

const char* alphabet = "ATGC";
inline char baseToChar (Base base) {
  return alphabet[base & 3];
}

inline Base getBase (Kmer kmer, Pos pos) {
  return (kmer >> ((pos - 1) << 1)) & 3;
}

inline Kmer setBase (Kmer kmer, Pos pos, Base base) {
  const int shift = (pos - 1) << 1;
  return (kmer & (((Kmer) -1) ^ (3 << shift))) | (base << shift);
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
  Kmer prefix = (kmer << 2) & kmerMask(len);
  iota (outgoing.begin(), outgoing.end(), prefix);
}

inline void getIncoming (Kmer kmer, Pos len, vguard<Kmer>& incoming) {
  Assert (incoming.size() == 4, "oops");
  Kmer prefix = kmer >> 2;
  const int shift = (len - 1) << 1;
  for (Base b = 0; b < 4; ++b)
    incoming[b] = prefix | (b << shift);
}

inline KmerFlags outgoingEdgeFlags (Kmer kmer, Pos len, vguard<Kmer>& outgoing, const vguard<KmerFlags>& flags) {
  getOutgoing (kmer, len, outgoing);
  KmerFlags f = 0;
  for (size_t n = 0; n < 4; ++n)
    if (flags[outgoing[n]] & KmerValid)
      f = f | (1 << n);
  return f;
}

int flagsToCountLookup[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
inline int flagsToCount (KmerFlags flags) {
  return flagsToCountLookup[flags & 0xf];
}

inline void pruneDeadEnds (Kmer kmer, Pos len, vguard<Kmer>& outgoing, vguard<KmerFlags>& flags) {
  if (flags[kmer] & KmerValid) {
    const KmerFlags outgoingFlags = outgoingEdgeFlags(kmer,len,outgoing,flags);
    if (outgoingFlags == 0) {
      LogThisAt(9,"Pruning " << kmerString(kmer,len) << endl);
      flags[kmer] = 0;
      vguard<Kmer> incoming (4);
      getIncoming (kmer, len, incoming);
      for (auto kmerIn: incoming)
      	pruneDeadEnds (kmerIn, len, outgoing, flags);
    } else
      LogThisAt(9,"Keeping " << kmerString(kmer,len) << " with " << flagsToCount(outgoingFlags) << " outgoing edges" << endl);
  }
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

struct MutPosIterator {
  const Pos len;
  const int maxDistance;
  vguard<Pos> pos;
  int nPos;
  vguard<Kmer> neighbors;
  int transversionDistance;
  MutPosIterator (Pos len, int maxDistance)
    : len(len),
      maxDistance(maxDistance),
      pos(maxDistance+1),
      transversionDistance(2)
  { }
  inline void reset() {
    pos[0] = 1;
    nPos = 1;
  }
  inline bool done() const {
    return nPos > maxDistance;
  }
  inline void next() {
    int nInc = nPos - 1;
    while (nInc >= 0)
      if (++pos[nInc] + (nPos - nInc - 1) <= len) {
	iota (pos.begin() + nInc + 1, pos.begin() + nPos, pos[nInc] + 1);
	break;
      } else
	--nInc;
    if (nInc < 0) {
      ++nPos;
      iota (pos.begin(), pos.begin() + nPos, 1);
    }
  }
  void getNeighbors (Kmer seq) {
    neighbors.clear();
    int nNbrs = 1;
    for (int n = 0; n < nPos; ++n)
      nNbrs *= 3;
    for (int nbrIdx = 0; nbrIdx < nNbrs; ++nbrIdx) {
      int d = 0;
      int tmp = nbrIdx;
      Kmer nbr = seq;
      for (int n = 0; n < nPos; ++n) {
	const int delta = (tmp % 3) + 1;
	tmp = tmp / 3;
	d += (delta == 2) ? 1 : transversionDistance;
	nbr = setBase (nbr, pos[n], (getBase (seq, pos[n]) + delta) & 3);
      }
      if (d <= maxDistance)
	neighbors.push_back (nbr);
    }
  }
};

int main (int argc, char** argv) {

  try {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "display this help message")
      ("kmerlen,k", po::value<int>()->default_value(15), "length of k-mers in de Bruijn graph")
      ("tandem,t", po::value<int>()->default_value(3), "reject local tandem duplications up to this length")
      ("invrep,i", po::value<int>()->default_value(5), "reject inverted repeats of this length (separated by at least 2 bases)")
      ("nbrdist,d", po::value<int>()->default_value(3), "exclude neighbors at up to this (transition+2*transversion)-distance")
      ("verbose,v", po::value<int>()->default_value(1), "verbosity level")
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
    const int nbrDist = vm.at("nbrdist").as<int>();

    const Kmer maxKmer = kmerMask(len);

    // Remove repeats
    ProgressLog (plogReps, 1);
    plogReps.initProgress ("Filtering %d-mer repeats", len);
    vguard<KmerFlags> kmerFlags (maxKmer+1);
    unsigned long long nKmersWithoutReps = 0;
    list<Kmer> kmersWithoutReps;
    for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
      plogReps.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      
      if (!hasExactTandemRepeat(kmer,len,maxTandemRepeatLen)
	  && !hasExactInvertedRepeat(kmer,len,invertedRepeatLen)) {
	LogThisAt(9,"Accepting " << kmerString(kmer,len) << endl);
	kmerFlags[kmer] = KmerValid;
	kmersWithoutReps.push_back (kmer);
	++nKmersWithoutReps;
      }
    }
    LogThisAt(2,"Found " << nKmersWithoutReps << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutReps/(1.+(double)maxKmer) << "%)" << endl);

    // Remove neighbors
    ProgressLog (plogNbrs, 1);
    plogNbrs.initProgress ("Rejecting neighbors up to (transition+2*transversion)-distance %d", nbrDist);
    unsigned long long nbrScans = 0, nbrsRejected = 0, nKmersWithoutNeighbors = 0;
    MutPosIterator mutPosIter (len, nbrDist);
    list<Kmer> kmersWithoutNeighbors;
    for (auto kmer: kmersWithoutReps) {
      ++nbrScans;
      plogNbrs.logProgress (nbrScans / (double) nKmersWithoutReps, "sequence %llu/%llu", nbrScans, nKmersWithoutReps);
      if (kmerFlags[kmer] & KmerValid) {
	kmersWithoutNeighbors.push_back (kmer);
	++nKmersWithoutNeighbors;
	for (mutPosIter.reset(); !mutPosIter.done(); mutPosIter.next()) {
	  mutPosIter.getNeighbors (kmer);
	  for (auto nbr: mutPosIter.neighbors)
	    if (kmerFlags[nbr] & KmerValid) {
	      LogThisAt(5,"Rejecting " << kmerString(kmer,len) << " neighbor " << kmerString(nbr,len) << endl);
	      kmerFlags[nbr] = 0;
	    }
	}
      }
    }
    LogThisAt(2,"Neighbor elimination left " << nKmersWithoutNeighbors << " candidate " << len << "-mers (" << setprecision(2) << 100*(double)nKmersWithoutNeighbors/(1.+(double)maxKmer) << "%)" << endl);

    // Remove dead ends
    ProgressLog (plogPrune, 1);
    plogNbrs.initProgress ("Pruning dead ends");
    vguard<Kmer> outgoing (4);
    for (auto kmer: kmersWithoutNeighbors)
      pruneDeadEnds (kmer, len, outgoing, kmerFlags);
    unsigned long long nPruned = 0, nUnpruned = 0;
    list<Kmer> unprunedKmers;
    for (auto kmer: kmersWithoutNeighbors) {
      ++nPruned;
      plogNbrs.logProgress (nPruned / (double) nKmersWithoutNeighbors, "sequence %llu/%llu", nPruned, nKmersWithoutNeighbors);
      if (kmerFlags[kmer] & KmerValid) {
	unprunedKmers.push_back (kmer);
	++nUnpruned;
      }
    }
    LogThisAt(2,"Dead-end pruning left " << nUnpruned << " candidate " << len << "-mers (" << setprecision(2) << 100*(double)nUnpruned/(1.+(double)maxKmer) << "%)" << endl);
    
  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}
