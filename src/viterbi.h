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
  vguard<Base> leftContext;
  vguard<IncomingTransScore> emit, null;
};

struct MachineScores {
  vguard<StateScores> stateScores;
  MachineScores (const Machine& machine, const InputModel& inputModel);
};

struct MutatorParams {
  double pDelOpen, pDelExtend, pTanDup, pTransition, pTransversion;
  Pos maxLen;
  bool local;
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
  size_t maxDupLen;
  vguard<LogProb> cell;

  static inline size_t nCells (size_t seqLen) {
    return (maxDupLen + 2) * (seqLen + 1);
  };

  inline size_t sCellIndex (Pos pos) const {
    return (maxDupLen + 2) * pos;
  };
  inline size_t dCellIndex (Pos pos) const {
    return (maxDupLen + 2) * pos + 1;
  };
  inline size_t fCellIndex (Pos pos, Pos idx) const {
    return (maxDupLen + 2) * pos + 1 + idx;
  };

public:
  const Machine& machine;
  const InputModel& inputModel;
  const MutatorParams& mutatorParams;
  const TokSeq& seq;
  MachineScores machineScores;
  MutatorScores mutatorScores;

  MutatorMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const TokSeq& seq);

  inline LogProb& sCell (Pos pos) { return cell[sCellIndex(pos)]; }
  inline LogProb& dCell (Pos pos) { return cell[dCellIndex(pos)]; }
  inline LogProb& fCell (Pos pos, Pos idx) { return cell[fCellIndex(pos,idx)]; }
  inline const LogProb sCell const (Pos pos) { return cell[sCellIndex(pos)]; }
  inline const LogProb dCell const (Pos pos) { return cell[dCellIndex(pos)]; }
  inline const LogProb fCell const (Pos pos, Pos idx) { return cell[fCellIndex(pos,idx)]; }
};

struct ViterbiMatrix : MutatorMatrix {

  ViterbiMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const TokSeq& seq);

  string traceback() const;
};

#endif /* VITERBI_INCLUDED */
