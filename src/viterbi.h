#ifndef VITERBI_INCLUDED
#define VITERBI_INCLUDED

#include "mutator.h"

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

class ViterbiMatrix {
private:
  size_t maxDupLen, nStates, seqLen;
  vguard<LogProb> cell;

  static inline size_t nCells (const Machine& machine, const TokSeq& seq) {
    return (maxDupLen + 2) * machine.nStates() * (seq.size() + 1);
  };

  inline size_t sCellIndex (State state, Pos pos) const {
    return (maxDupLen + 2) * (pos * nStates + state);
  };
  inline size_t dCellIndex (State state, Pos pos) const {
    return (maxDupLen + 2) * (pos * nStates + state) + 1;
  };
  inline size_t tCellIndex (State state, Pos pos, Pos idx) const {
    return (maxDupLen + 2) * (pos * nStates + state) + 1 + idx;
  };

  inline LogProb& sCell (State state, Pos pos) { return cell[sCellIndex(state,pos)]; }
  inline LogProb& dCell (State state, Pos pos) { return cell[dCellIndex(state,pos)]; }
  inline LogProb& tCell (State state, Pos pos, Pos idx) { return cell[tCellIndex(state,pos,idx)]; }

public:
  const Machine& machine;
  const InputModel& inputModel;
  const MutatorParams& mutatorParams;
  const TokSeq& seq;
  const MachineScores machineScores;
  const MutatorScores mutatorScores;

  ViterbiMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const TokSeq& seq);
  string traceback() const;

  inline const LogProb sCell const (State state, Pos pos) { return cell[sCellIndex(state,pos)]; }
  inline const LogProb dCell const (State state, Pos pos) { return cell[dCellIndex(state,pos)]; }
  inline const LogProb tCell const (State state, Pos pos, Pos idx) { return cell[tCellIndex(state,pos,idx)]; }
};

#endif /* VITERBI_INCLUDED */
