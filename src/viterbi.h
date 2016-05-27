#ifndef VITERBI_INCLUDED
#define VITERBI_INCLUDED

#include "mutator.h"
#include "fastseq.h"

struct InputModel {
  map<InputSymbol,double> symProb;
  InputModel (const string& inputAlphabet, double controlProb);
  InputModel();
};

struct IncomingTransScore {
  State src;
  LogProb score;
  InputSymbol in;
  Base base;
};

struct StateScores {
  vguard<Base> leftContext;
  vguard<IncomingTransScore> emit, null;
  inline Base base() const { return leftContext.back(); }
};

struct MachineScores {
  vguard<StateScores> stateScores;
  MachineScores (const Machine& machine, const InputModel& inputModel);
};

class ViterbiMatrix {
private:
  size_t maxDupLen, nStates, seqLen;
  vguard<LogProb> cell;

  static inline size_t nCells (const Machine& machine, const MutatorParams& params, const TokSeq& seq) {
    return (params.maxDupLen() + 2) * machine.nStates() * (seq.size() + 1);
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
  const FastSeq& fastSeq;
  const TokSeq seq;
  const MachineScores machineScores;
  const MutatorScores mutatorScores;

  LogProb loglike;
  
  ViterbiMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const FastSeq& fastSeq);
  string traceback() const;

  inline LogProb sCell (State state, Pos pos) const { return cell[sCellIndex(state,pos)]; }
  inline LogProb dCell (State state, Pos pos) const { return cell[dCellIndex(state,pos)]; }
  inline LogProb tCell (State state, Pos pos, Pos idx) const { return cell[tCellIndex(state,pos,idx)]; }

  inline Pos maxDupLenAt (const StateScores& ss) const { return min ((Pos) maxDupLen, (Pos) ss.leftContext.size()); }
  inline Base tanDupBase (const StateScores& ss, Pos dupIdx) const { return ss.leftContext[ss.leftContext.size() - 1 - dupIdx]; }
};

#endif /* VITERBI_INCLUDED */
