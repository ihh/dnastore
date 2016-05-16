#ifndef TRANSBUILDER_INCLUDED
#define TRANSBUILDER_INCLUDED

#include <string>
#include "vguard.h"
#include "util.h"
#include "kmer.h"
#include "pattern.h"

using namespace std;

struct TransBuilder {
  static vguard<int> edgeFlagsToCountLookup;

  // specified at creation
  const Pos len;
  const Kmer maxKmer;

  // config
  Pos maxTandemRepeatLen, invertedRepeatLen;
  set<KmerLen> excludedMotif, excludedMotifRevComp;
  set<KmerLen> sourceMotif;
  int controlWords;
  
  // work variables
  vguard<bool> kmerValid;
  list<Kmer> kmers;
  map<Kmer,State> kmerState;
  map<Kmer,int> timesAvoided;
  
  TransBuilder (Pos len);

  void removeRepeats();
  void pruneUnreachable();
  void pruneDeadEnds();
  void indexStates();
  
  void output (ostream& out);

  void doDFS (Kmer kmer, map<Kmer,int>& distance, bool forwards) const;

  set<Kmer> neighbors (const set<Kmer>& start, int steps, bool forwards) const;
  set<Kmer> neighbors (Kmer start, int steps, bool forwards) const;

  int stepsToReach (KmerLen motif, int maxSteps = 64) const;

  inline void pruneDeadEnds (Kmer kmer) {
    EdgeVector in, out;
    if (kmerValid[kmer] && !endsWithMotif(kmer,len,sourceMotif)) {
      const EdgeFlags inFlags = incomingEdgeFlags(kmer,in);
      const EdgeFlags outFlags = outgoingEdgeFlags(kmer,out);
      const bool prune = inFlags == 0 || outFlags == 0;
      LogThisAt(9,(prune ? "Pruning" : "Keeping") << " " << kmerString(kmer,len) << " with " << edgeFlagsToCount(inFlags) << " incoming and " << edgeFlagsToCount(outFlags) << " outgoing edges" << endl);
      if (prune) {
	kmerValid[kmer] = 0;
	for (auto kmerIn: in)
	  pruneDeadEnds (kmerIn);
	for (auto kmerOut: out)
	  pruneDeadEnds (kmerOut);
      }
    }
  }
 
  inline void getOutgoing (Kmer kmer, EdgeVector& outgoing) const {
    Assert (outgoing.size() == 4, "oops");
    const Kmer prefix = (kmer << 2) & kmerMask(len);
    iota (outgoing.begin(), outgoing.end(), prefix);
  }

  inline void getIncoming (Kmer kmer, EdgeVector& incoming) const {
    Assert (incoming.size() == 4, "oops");
    const Kmer prefix = kmer >> 2;
    const int shift = (len - 1) << 1;
    for (Base b = 0; b < 4; ++b)
      incoming[b] = prefix | (((Kmer) b) << shift);
  }

  inline EdgeFlags outgoingEdgeFlags (Kmer kmer, EdgeVector& outgoing) {
    getOutgoing (kmer, outgoing);
    EdgeFlags f = 0;
    for (size_t n = 0; n < 4; ++n)
      if (kmerValid[outgoing[n]] && !endsWithMotif(outgoing[n],len,sourceMotif))
	f = f | (1 << n);
    return f;
  }

  inline EdgeFlags incomingEdgeFlags (Kmer kmer, EdgeVector& incoming) {
    getIncoming (kmer, incoming);
    EdgeFlags f = 0;
    for (size_t n = 0; n < 4; ++n)
      if (kmerValid[incoming[n]])
	f = f | (1 << n);
    return f;
  }

  inline int edgeFlagsToCount (EdgeFlags flags) {
    return edgeFlagsToCountLookup[flags & 0xf];
  }

  inline bool betterDest (Kmer x, Kmer y) {  // true if x is preferred target
    const int xa = timesAvoided[x], ya = timesAvoided[y];
    const double xgc = gcNonuniformity(x,len);
    const double ygc = gcNonuniformity(y,len);
    return xa == ya
      ? (xgc == ygc
	 ? (kmerEntropy(x,len) >= kmerEntropy(y,len))
	 : (xgc < ygc))
      : (xa < ya);
  }
};

#endif /* TRANSBUILDER_INCLUDED */
