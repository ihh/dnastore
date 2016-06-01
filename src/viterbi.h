#ifndef VITERBI_INCLUDED
#define VITERBI_INCLUDED

#include "mutator.h"
#include "fastseq.h"

struct InputModel {
  string inputAlphabet;
  map<InputSymbol,double> symProb;
  InputModel (const string& inputAlphabet, double symWeight = 1., double controlWeight = 1.);
  string toString() const;
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
  typedef size_t MutStateIndex;
  size_t maxDupLen, nStates, seqLen;
  vguard<LogProb> cell;

  static inline size_t nCells (const Machine& machine, const MutatorParams& params, const FastSeq& seq) {
    return (params.maxDupLen() + 2) * machine.nStates() * (seq.length() + 1);
  };

  inline MutStateIndex sMutStateIndex() const { return 0; }
  inline MutStateIndex dMutStateIndex() const { return 1; }
  inline MutStateIndex tMutStateIndex (Pos dupIdx) const { return 2 + dupIdx; }

  inline Pos tMutStateDupIdx (MutStateIndex mutStateIndex) const { return mutStateIndex - 2; }
  inline bool isTMutStateIndex (MutStateIndex mutStateIndex) const {
    return mutStateIndex >= tMutStateIndex(0) && mutStateIndex <= tMutStateIndex(maxDupLen-1);
  }

  inline string mutStateName (MutStateIndex mutStateIndex) const {
    return mutStateIndex==0 ? string("S") : (mutStateIndex==1 ? string("D") : (string("T") + to_string(tMutStateDupIdx(mutStateIndex)+1)));
  }
  
  inline size_t cellIndex (State state, Pos pos, MutStateIndex mutState) const {
    return (maxDupLen + 2) * (pos * nStates + state) + mutState;
  };
  inline size_t sCellIndex (State state, Pos pos) const {
    return cellIndex (state, pos, sMutStateIndex());
  };
  inline size_t dCellIndex (State state, Pos pos) const {
    return cellIndex (state, pos, dMutStateIndex());
  };
  inline size_t tCellIndex (State state, Pos pos, Pos dupIdx) const {
    return cellIndex (state, pos, tMutStateIndex(dupIdx));
  };

  inline LogProb& sCell (State state, Pos pos) { return cell[sCellIndex(state,pos)]; }
  inline LogProb& dCell (State state, Pos pos) { return cell[dCellIndex(state,pos)]; }
  inline LogProb& tCell (State state, Pos pos, Pos idx) { return cell[tCellIndex(state,pos,idx)]; }

  inline LogProb getCell (State state, Pos pos, MutStateIndex mutState) const { return cell[cellIndex(state,pos,mutState)]; }
  inline LogProb& loglike() { return sCell (machine.nStates() - 1, seqLen); }
  
public:
  const Machine& machine;
  const InputModel& inputModel;
  const MutatorParams& mutatorParams;
  const FastSeq& fastSeq;
  const TokSeq seq;
  const MachineScores machineScores;
  const MutatorScores mutatorScores;

  ViterbiMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const FastSeq& fastSeq);
  string toString() const;
  string traceback() const;

  inline LogProb sCell (State state, Pos pos) const { return cell[sCellIndex(state,pos)]; }
  inline LogProb dCell (State state, Pos pos) const { return cell[dCellIndex(state,pos)]; }
  inline LogProb tCell (State state, Pos pos, Pos dupIdx) const { return cell[tCellIndex(state,pos,dupIdx)]; }

  inline LogProb loglike() const { return sCell (machine.nStates() - 1, seqLen); }
  
  inline Pos maxDupLenAt (const StateScores& ss) const { return min ((Pos) maxDupLen, (Pos) ss.leftContext.size()); }
  inline Base tanDupBase (const StateScores& ss, Pos dupIdx) const { return ss.leftContext[ss.leftContext.size() - 1 - dupIdx]; }
};

vguard<FastSeq> decodeFastSeqs (const char* filename, const Machine& machine, const MutatorParams& mutatorParams);

#endif /* VITERBI_INCLUDED */
