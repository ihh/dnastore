#include <list>
#include "viterbi.h"

InputModel::InputModel (const string& inputAlphabet, double controlProb) {
  size_t nControls = 0;
  for (char c: inputAlphabet)
    if (Machine::isControl(c))
      ++nControls;
  for (char c: inputAlphabet)
    symProb[c] = Machine::isControl(c)
      ? (controlProb / (double) nControls)
      : ((1. - controlProb) / (double) (inputAlphabet.size() - nControls));
}

MachineScores::MachineScores (const Machine& machine, const InputModel& inputModel)
  : stateScores (machine.nStates())
{
  machine.verifyContexts();
  Assert (machine.isWaitingMachine(), "Not a waiting machine");
  for (char c: machine.outputAlphabet())
    Assert (isValidToken(c,dnaAlphabetString), "Not a DNA-outputting machine");

  for (State s = 0; s < machine.nStates(); ++s) {
    const MachineState& ms = machine.state[s];
    StateScores& ss = stateScores[s];
    ss.leftContext.reserve (ms.leftContext.size());
    for (size_t n = 0; n < ms.leftContext.size(); ++n)
      ss.leftContext.push_back (charToBase (ms.leftContext[n]));
    for (const auto& t: ms.trans) {
      if (t.inputEmpty() || t.isEOF() || inputModel.symProb.count(t.in)) {
	IncomingTransScore its;
	its.src = s;
	its.score = inputModel.symProb.count(t.in) ? log(inputModel.symProb.at(t.in)) : 0;
	its.in = t.in;
	StateScores& destStateScores = stateScores[t.dest];
	if (t.outputEmpty())
	  destStateScores.null.push_back (its);
	else {
	  its.base = charToBase (t.out);
	  destStateScores.emit.push_back (its);
	}
      }
    }
  }
}

ViterbiMatrix::ViterbiMatrix (const Machine& machine, const InputModel& inputModel, const MutatorParams& mutatorParams, const FastSeq& fastSeq)
  : maxDupLen (mutatorParams.maxDupLen()),
    nStates (machine.nStates()),
    seqLen (fastSeq.length()),
    cell (nCells (machine, mutatorParams, seq), -numeric_limits<double>::infinity()),
    machine (machine),
    inputModel (inputModel),
    mutatorParams (mutatorParams),
    fastSeq (fastSeq),
    seq (fastSeq.tokens (dnaAlphabetString)),
    machineScores (machine, inputModel),
    mutatorScores (mutatorParams),
    loglike (-numeric_limits<double>::infinity())
{
  for (State state = 0; state < machine.nStates(); ++state) {
    sCell(state,0) = 0;
    for (Pos pos = 1; pos <= seqLen; ++pos) {
      const Base base = seq[pos-1];
      const StateScores& ss = machineScores.stateScores[state];
      const auto mdl = maxDupLenAt(ss);

      for (const auto& its: ss.emit) {
	dCell(state,pos) = max (dCell(state,pos),
				max (dCell(its.src,pos) + its.score + mutatorScores.delExtend,
				     sCell(its.src,pos) + its.score + mutatorScores.delOpen));
	sCell(state,pos) = max (sCell(state,pos),
				sCell(its.src,pos-1) + its.score + mutatorScores.noGap + mutatorScores.sub[its.base][base]);
      }
      for (const auto& its: ss.null) {
	dCell(state,pos) = max (dCell(state,pos),
				dCell(its.src,pos) + its.score);
	sCell(state,pos) = max (sCell(state,pos),
				sCell(its.src,pos) + its.score);
      }
      sCell(state,pos) = max (sCell(state,pos),
			      dCell(state,pos) + mutatorScores.delEnd);

      if (mdl > 0) {
	sCell(state,pos) = max (sCell(state,pos),
				tCell(state,pos-1,0) + mutatorScores.sub[tanDupBase(ss,0)][base]);

	for (Pos dupIdx = 0; dupIdx < mdl - 1; ++dupIdx)
	  tCell(state,pos,dupIdx) = tCell(state,pos-1,dupIdx+1) + mutatorScores.sub[tanDupBase(ss,dupIdx+1)][base];

	for (Pos dupIdx = 0; dupIdx < mdl; ++dupIdx)
	  tCell(state,pos,dupIdx) = sCell(state,pos) + mutatorScores.tanDup + mutatorScores.len[dupIdx];
      }
    }

    loglike = max (loglike, sCell(state,seqLen));
  }
}

string ViterbiMatrix::traceback() const {
  list<char> trace;
  
  LogProb best;
  State bestState;
  Pos bestPos;
  InputSymbol bestInSym;
  MutStateIndex bestMutState;

  auto updateBest = [&] (State srcState, Pos srcPos, MutStateIndex srcMutState, LogProb transScore, InputSymbol inSym) -> void {
    const LogProb score = getCell(srcState,srcPos,srcMutState) + transScore;
    if (score > best) {
      best = score;
      bestState = srcState;
      bestPos = srcPos;
      bestMutState = srcMutState;
      bestInSym = inSym;
    }
  };

  auto checkBest = [&]() -> void {
    const LogProb expected = getCell(bestState,bestPos,bestMutState);
    Assert (abs((best - expected) / expected) < 1e-6, "Traceback failure");
  };

  best = -numeric_limits<double>::infinity();
  for (State s = 0; s < machine.nStates(); ++s)
    updateBest (s, seqLen, sMutStateIndex(), 0, MachineNull);
  checkBest();

  while (bestPos >= 0) {
    const Base base = seq[bestPos-1];
    const StateScores& ss = machineScores.stateScores[bestState];
    const auto mdl = maxDupLenAt(ss);
    if (bestMutState == sMutStateIndex()) {

      for (const auto& its: ss.emit)
	updateBest (its.src, bestPos-1, sMutStateIndex(), its.score + mutatorScores.noGap + mutatorScores.sub[its.base][base], its.in);
      for (const auto& its: ss.null)
	updateBest (its.src, bestPos, sMutStateIndex(), its.score, MachineNull);
      updateBest (bestState, bestPos, dMutStateIndex(), mutatorScores.delEnd, MachineNull);

      if (mdl > 0)
	updateBest (bestState, bestPos-1, tMutStateIndex(0), mutatorScores.sub[tanDupBase(ss,0)][base], MachineNull);

    } else if (bestMutState == dMutStateIndex()) {

      for (const auto& its: ss.emit) {
	updateBest (its.src, bestPos, dMutStateIndex(), its.score + mutatorScores.delExtend, its.in);
	updateBest (its.src, bestPos, sMutStateIndex(), its.score + mutatorScores.delOpen, its.in);
      }
      for (const auto& its: ss.null)
	updateBest (its.src, bestPos, dMutStateIndex(), its.score, MachineNull);

    } else if (isTMutStateIndex(bestMutState)) {

      const Pos dupIdx = tMutStateDupIdx (bestMutState);
      if (dupIdx < mdl - 1)
	updateBest (bestState, bestPos-1, dupIdx+1, mutatorScores.sub[tanDupBase(ss,dupIdx+1)][base], MachineNull);
      updateBest (bestState, bestPos, sMutStateIndex(), mutatorScores.tanDup + mutatorScores.len[dupIdx], MachineNull);

    } else
      Abort ("Unknown traceback state");

    checkBest();
    if (bestInSym)
      trace.push_front (bestInSym);
  }

  return string (trace.begin(), trace.end());
}
