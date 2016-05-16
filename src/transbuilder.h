#ifndef TRANSBUILDER_INCLUDED
#define TRANSBUILDER_INCLUDED

#include <string>
#include "vguard.h"
#include "util.h"
#include "kmer.h"
#include "pattern.h"

using namespace std;

struct TransBuilder {
  const Pos len;
  const Kmer maxKmer;
  vguard<bool> kmerValid;
  Pos maxTandemRepeatLen, invertedRepeatLen;
  set<KmerLen> excludedMotif, excludedMotifRevComp;
  set<KmerLen> sourceMotif;

  list<Kmer> kmers;
  map<Kmer,State> kmerState;
  static vguard<int> edgeFlagsToCountLookup;

  TransBuilder (Pos len);

  void removeRepeats();
  void keepDFS();
  void pruneDeadEnds();
  void indexStates();
  
  void output (ostream& out);

  void doDFS (Kmer kmer, set<Kmer>& seen) const;

  inline void pruneDeadEnds (Kmer kmer) {
    vguard<Kmer> in (4), out (4);
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
 
  inline void getOutgoing (Kmer kmer, vguard<Kmer>& outgoing) const {
    Assert (outgoing.size() == 4, "oops");
    const Kmer prefix = (kmer << 2) & kmerMask(len);
    iota (outgoing.begin(), outgoing.end(), prefix);
  }

  inline void getIncoming (Kmer kmer, vguard<Kmer>& incoming) const {
    Assert (incoming.size() == 4, "oops");
    const Kmer prefix = kmer >> 2;
    const int shift = (len - 1) << 1;
    for (Base b = 0; b < 4; ++b)
      incoming[b] = prefix | (((Kmer) b) << shift);
  }

  inline EdgeFlags outgoingEdgeFlags (Kmer kmer, vguard<Kmer>& outgoing) {
    getOutgoing (kmer, outgoing);
    EdgeFlags f = 0;
    for (size_t n = 0; n < 4; ++n)
      if (kmerValid[outgoing[n]] && !endsWithMotif(outgoing[n],len,sourceMotif))
	f = f | (1 << n);
    return f;
  }

  inline EdgeFlags incomingEdgeFlags (Kmer kmer, vguard<Kmer>& incoming) {
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

};

#endif /* TRANSBUILDER_INCLUDED */
