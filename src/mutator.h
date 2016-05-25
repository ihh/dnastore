#ifndef MUTATOR_INCLUDED
#define MUTATOR_INCLUDED

#include <iostream>
#include "kmer.h"
#include "trans.h"
#include "logsumexp.h"

struct MutatorParams {
  double pDelOpen, pDelExtend, pTanDup, pTransition, pTransversion;
  vguard<double> pLen;
  bool local;

  void writeJSON (ostream& out) const;
  void readJSON (istream& in);
  static MutatorParams fromJSON (istream& in);
  static MutatorParams fromFile (const char* filename);

  inline double pMatch() const { return 1. - pTransition - pTransversion; }
  inline double pSub (Base x, Base y) const {
    return x == y ? pMatch() : (isTransition(x,y) ? pTransition : (pTransversion/2.));
  }
  inline double pFwdDup() const { return 0.; }
  inline double pRevDup() const { return 0.; }
  inline double pNoGap() const { return 1. - pDelOpen - pTanDup; }
  inline double pDelEnd() const { return 1. - pDelExtend; }
  inline size_t maxLen() const { return pLen.size(); }
};

struct MutatorScores {
  LogProb delOpen, tanDup, noGap;
  LogProb delExtend, delEnd;
  vguard<vguard<LogProb> > sub;  // sub[base][observed]
  vguard<LogProb> len;
  MutatorScores (const MutatorParams& params);
};

struct MutatorCounts {
  double nDelOpen, nTanDup, nNoGap;
  vguard<vguard<double> > nSub;
  vguard<double> nLen;

  MutatorCounts& initMaxLen (size_t maxLen);
  MutatorCounts& initLaplace (double n = 1);

  MutatorParams mlParams() const;

  LogProb logPrior (const MutatorParams& params) const;
  LogProb logLikelihood (const MutatorParams& params) const;
};

#endif /* MUTATOR_INCLUDED */
