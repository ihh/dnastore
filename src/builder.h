#ifndef BUILDER_INCLUDED
#define BUILDER_INCLUDED

#include <string>
#include "vguard.h"
#include "util.h"
#include "kmer.h"
#include "pattern.h"
#include "trans.h"

using namespace std;

typedef unsigned char EdgeFlags;

struct TransBuilder {
  static vguard<int> edgeFlagsToCountLookup;

  // specified at creation
  const Pos len;
  const Kmer maxKmer;

  // config
  Pos maxTandemRepeatLen, invertedRepeatLen;
  set<KmerLen> excludedMotif, excludedMotifRevComp;
  set<KmerLen> sourceMotif;
  bool keepDegenerates;
  size_t nControlWords;
  
  // work variables
  vguard<bool> kmerValid;
  list<Kmer> kmers;
  vguard<Kmer> controlWord;
  vguard<Pos> controlWordSteps;
  vguard<map<Kmer,list<Kmer> > > controlWordPath;
  vguard<vguard<set<Kmer> > > controlWordIntermediates;
  map<Kmer,EdgeFlags> kmerOutFlags;
  set<pair<Kmer,Kmer> > droppedEdge;

  State nStates, nCodingStates;
  map<Kmer,State> kmerState, kmerStateZero, kmerStateOne;
  vguard<vguard<map<Kmer,State> > > controlKmerState;
  
  TransBuilder (Pos len);

  void findCandidates();
  void pruneUnreachable();
  void pruneDeadEnds();
  void buildEdges();
  void indexStates();

  Machine makeMachine();
  
  void assertKmersCorrect() const;
  
  void doDFS (Kmer kmer, map<Kmer,Pos>& distance) const;

  set<Kmer> kmersEndingWith (KmerLen motif) const;
  Pos stepsToReach (KmerLen motif, int maxSteps = 64) const;
  void getControlWords();
  
  map<Kmer,list<Kmer> > pathsTo (Kmer dest, int steps) const;
  MachineTransition controlTrans (State srcState, Kmer destKmer, size_t nControlWord, size_t step) const;
  Kmer nextIntermediateKmer (Kmer srcKmer, size_t nControlWord, size_t step) const;
  char controlChar (size_t nControlWord) const;
  
  inline void pruneDeadEnds (Kmer kmer) {
    EdgeVector in, out;
    if (kmerValid[kmer] && !endsWithMotif(kmer,len,sourceMotif)) {
      const int inCount = countIncoming(kmer,in), outCount = countOutgoing(kmer,out);
      const bool prune = inCount == 0 || outCount == 0;
      LogThisAt(9,(prune ? "Pruning" : "Keeping") << " " << kmerString(kmer,len) << " with " << inCount << " incoming and " << outCount << " outgoing edges" << endl);
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
      if (kmerValid[outgoing[n]]
	  && !endsWithMotif(outgoing[n],len,sourceMotif)
	  && !droppedEdge.count (pair<Kmer,Kmer> (kmer, outgoing[n])))
	f = f | (1 << n);
    return f;
  }

  inline EdgeFlags incomingEdgeFlags (Kmer kmer, EdgeVector& incoming) {
    getIncoming (kmer, incoming);
    EdgeFlags f = 0;
    for (size_t n = 0; n < 4; ++n)
      if (kmerValid[incoming[n]]
	  && !droppedEdge.count (pair<Kmer,Kmer> (incoming[n], kmer)))
	f = f | (1 << n);
    return f;
  }

  inline int edgeFlagsToCount (EdgeFlags flags) {
    return edgeFlagsToCountLookup[flags & 0xf];
  }

  inline int countOutgoing (Kmer kmer, EdgeVector& out) {
    return edgeFlagsToCount (outgoingEdgeFlags (kmer, out));
  }

  inline int countIncoming (Kmer kmer, EdgeVector& in) {
    return edgeFlagsToCount (incomingEdgeFlags (kmer, in));
  }

  inline int countIncoming (Kmer kmer) {
    EdgeVector e;
    return countIncoming (kmer, e);
  }

  inline int countOutgoing (Kmer kmer) {
    EdgeVector e;
    return countOutgoing (kmer, e);
  }

  inline bool betterDest (Kmer x, Kmer y) {  // true if x is preferred over y as destination state
    const int xi = countIncoming(x), yi = countIncoming(y);
    const double xgc = gcNonuniformity(x,len), ygc = gcNonuniformity(y,len);
    return xi == yi
      ? (xgc == ygc
	 ? (kmerEntropy(x,len) >= kmerEntropy(y,len))
	 : (xgc < ygc))
      : (xi < yi);
  }

  inline EdgeFlags dropWorseEdge (Kmer src, EdgeFlags flags, const EdgeVector& out, size_t edge1, size_t edge2) {
    const size_t e = betterDest(out[edge1],out[edge2]) ? edge2 : edge1;
    LogThisAt(4,"Dropping "
	      << (countIncoming(out[e]) == 1 ? "last " : "")
	      << "edge to " << kmerString(out[e],len)
	      << " from " << kmerString(src,len)
	      << endl);
    droppedEdge.insert (pair<Kmer,Kmer> (src, out[e]));
    return flags & (0xf ^ (1 << e));
  }
};

#endif /* BUILDER_INCLUDED */
