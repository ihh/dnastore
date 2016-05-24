#ifndef VITERBI_INCLUDED
#define VITERBI_INCLUDED

#include "trans.h"

struct InputModel {
  map<InputSymbol,double> symProb;
  InputModel (const string& inputAlphabet, double controlProb);
};

struct IncomingTransScore {
  State src;
  LogProb score;
};

struct StateScores {
  Base base;
  vguard<IncomingTransScore> emit, null;
};

struct MachineScores {
  vguard<StateScores> stateScores;
  MachineScores (const Machine& machine, const InputModel& inputModel);
};

struct MutatorParams {
  double pDelOpen, pDelExtend, pTanDup, pTransition, pTransversion;
  Pos maxLen;
  inline double pMatch() const { return 1. - pTransition - pTransversion; }
  inline double pSub (Base x, Base y) const {
    return x == y ? pMatch() : (isTransition(x,y) ? pTransition : (pTransversion/2.));
  }
  inline double pFwdDup() const { return 0.; }
  inline double pRevDup() const { return 0.; }
  inline double pNoGap() const { return 1. - pDelOpen - pTanDup; }
  inline double pDelEnd() const { return 1. - pDelExtend; }
  inline double pLen (Pos len) const { return 1. / (double) maxLen; }
};

struct MutatorScores {
  LogProb match, delOpen, tanDup, noGap;
  LogProb delExtend, delEnd;
  vguard<vguard<LogProb> > sub;  // sub[base][observed]
  vguard<LogProb> len;
  MutatorScores (const MutatorParams& params);
};

class MutatorMatrix {
private:
  vguard<LogProb> cell;
public:
  const Machine& machine;
  const MutatorScores& scores;
  const TokSeq& seq;
  MutatorMatrix (const MutatorScores& scores, const TokSeq& seq);
};

struct ViterbiMatrix {
};

#endif /* VITERBI_INCLUDED */
