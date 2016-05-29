#include <list>
#include <iomanip>
#include "viterbi.h"
#include "logger.h"

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

InputModel::InputModel() {
  symProb[MachineSOF] = 1;
  symProb[MachineEOF] = 1;
  symProb[MachineBit0] = 1;
  symProb[MachineBit1] = 1;
  for (char c = MachineControlFirst; c <= MachineControlLast; ++c)
    symProb[c] = 1;
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
    for (char lc: ms.leftContext)
      if (lc != MachineWildContext)
	ss.leftContext.push_back (charToBase (lc));
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
  : maxDupLen (min (machine.maxLeftContext(), mutatorParams.maxDupLen())),
    nStates (machine.nStates()),
    seqLen (fastSeq.length()),
    cell (nCells (machine, mutatorParams, fastSeq), -numeric_limits<double>::infinity()),
    machine (machine),
    inputModel (inputModel),
    mutatorParams (mutatorParams),
    fastSeq (fastSeq),
    seq (fastSeq.tokens (dnaAlphabetString)),
    machineScores (machine, inputModel),
    mutatorScores (mutatorParams)
{
  if (mutatorParams.local)
    for (State state = 0; state < machine.nStates(); ++state)
      sCell(state,0) = 0;
  else
    sCell(0,0) = 0;

  ProgressLog (plog, 2);
  plog.initProgress ("Filling Viterbi matrix (%d*%d cells)", seqLen, machine.nStates());

  for (Pos pos = 0; pos <= seqLen; ++pos) {
    plog.logProgress (pos / (double) seqLen, "row %d/%d", pos, seqLen);
    for (State state = 0; state < machine.nStates(); ++state) {
      const StateScores& ss = machineScores.stateScores[state];
      const auto mdl = maxDupLenAt(ss);

      for (const auto& its: ss.emit) {
	dCell(state,pos) = max (dCell(state,pos),
				max (dCell(its.src,pos) + its.score + mutatorScores.delExtend,
				     sCell(its.src,pos) + its.score + mutatorScores.delOpen));
	if (pos > 0)
	  sCell(state,pos) = max (sCell(state,pos),
				  sCell(its.src,pos-1) + its.score + mutatorScores.noGap + mutatorScores.sub[its.base][seq[pos-1]]);
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
	if (pos > 0) {
	  sCell(state,pos) = max (sCell(state,pos),
				  tCell(state,pos-1,0) + mutatorScores.sub[tanDupBase(ss,0)][seq[pos-1]]);

	  for (Pos dupIdx = 0; dupIdx < mdl - 1; ++dupIdx)
	    tCell(state,pos,dupIdx) = tCell(state,pos-1,dupIdx+1) + mutatorScores.sub[tanDupBase(ss,dupIdx+1)][seq[pos-1]];
	}

	for (Pos dupIdx = 0; dupIdx < mdl; ++dupIdx)
	  tCell(state,pos,dupIdx) = max (tCell(state,pos,dupIdx),
					 sCell(state,pos) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
      }
    }
  }

  if (mutatorParams.local)
    for (State state = 0; state < machine.nStates(); ++state)
      loglike() = max (loglike(), sCell(state,seqLen));

  LogThisAt(10,"Viterbi matrix:\n" << toString());
}

string ViterbiMatrix::toString() const {
  ostringstream out;
  const size_t sw = machine.stateNameWidth();
  for (Pos pos = 0; pos <= seqLen; ++pos) {
    for (State state = 0; state < machine.nStates(); ++state) {
      out << setw(4) << pos << " "
	  << setw(sw) << machine.state[state].name << " "
	  << setw(10) << setprecision(6) << sCell(state,pos) << "(S) "
	  << setw(10) << setprecision(6) << dCell(state,pos) << "(D) ";
      for (Pos i = 0; i < maxDupLen; ++i)
	out << setw(10) << setprecision(6) << tCell(state,pos,i) << "(T" << i+1 << ") ";
      out << "\n";
    }
  }
  return out.str();
}

string ViterbiMatrix::traceback() const {
  list<char> trace;

  if (!(loglike() > -numeric_limits<double>::infinity())) {
    Warn ("No valid Viterbi decoding found");
    return "";
  }
  
  State state = machine.nStates() - 1, bestState;
  Pos pos = seqLen, bestPos;
  MutStateIndex mutState = 0, bestMutState;
  LogProb best;
  InputSymbol bestInSym;
  const IncomingTransScore* bestIts;
  bool foundBest;
  
  auto initBest = [&]() -> void {
    best = -numeric_limits<double>::infinity();
    foundBest = false;
    LogThisAt (9, "Traceback at (" << machine.state[state].name << "," << pos << "," << mutStateName(mutState) << ")" << endl);
  };

  auto updateBest = [&] (State srcState, Pos srcPos, MutStateIndex srcMutState, LogProb transScore, const IncomingTransScore* its) -> void {
    const LogProb score = getCell(srcState,srcPos,srcMutState) + transScore;
    if (score > best) {
      best = score;
      bestState = srcState;
      bestPos = srcPos;
      bestMutState = srcMutState;
      bestInSym = its ? its->in : MachineNull;
      bestIts = its;
      foundBest = true;
    }
  };

  auto checkBest = [&]() -> void {
    const LogProb expected = getCell(state,pos,mutState);
    Assert (abs((best - expected) / (abs(expected) < 1e-6 ? 1 : expected)) < 1e-6, "Traceback failure at (%s,%d,%s): computed traceback score (%g) didn't match stored value in matrix (%g)", machine.state[state].name.c_str(), pos, mutStateName(mutState).c_str(), best, expected);
    Assert (foundBest, "Traceback failure at (%s,%d,%s): couldn't find source state", machine.state[state].name.c_str(), pos, mutStateName(mutState).c_str());
    state = bestState;
    pos = bestPos;
    mutState = bestMutState;
  };

  initBest();
  if (mutatorParams.local)
    for (State s = 0; s < machine.nStates(); ++s)
      updateBest (s, seqLen, sMutStateIndex(), 0, NULL);
  else
    updateBest (machine.nStates() - 1, seqLen, sMutStateIndex(), 0, NULL);
  checkBest();

  while (pos >= 0 && state > 0) {
    const StateScores& ss = machineScores.stateScores[state];
    const auto mdl = maxDupLenAt(ss);
    initBest();
    if (mutState == sMutStateIndex()) {

      if (pos > 0)
	for (const auto& its: ss.emit)
	  updateBest (its.src, pos-1, sMutStateIndex(), its.score + mutatorScores.noGap + mutatorScores.sub[its.base][seq[pos-1]], &its);
      for (const auto& its: ss.null)
	updateBest (its.src, pos, sMutStateIndex(), its.score, &its);
      updateBest (state, pos, dMutStateIndex(), mutatorScores.delEnd, NULL);

      if (mdl > 0 && pos > 0)
	updateBest (state, pos-1, tMutStateIndex(0), mutatorScores.sub[tanDupBase(ss,0)][seq[pos-1]], NULL);

      if (pos == 0 && mutatorParams.local)
	updateBest (0, 0, sMutStateIndex(), 0, NULL);

      if (bestIts && bestPos < pos && seq[pos-1] != bestIts->base)
	LogThisAt(3,"Substitution at " << pos-1 << ": " << baseToChar(bestIts->base) << " -> " << baseToChar(seq[pos-1]) << endl);

    } else if (mutState == dMutStateIndex()) {

      for (const auto& its: ss.emit) {
	updateBest (its.src, pos, dMutStateIndex(), its.score + mutatorScores.delExtend, &its);
	updateBest (its.src, pos, sMutStateIndex(), its.score + mutatorScores.delOpen, &its);
      }
      for (const auto& its: ss.null)
	updateBest (its.src, pos, dMutStateIndex(), its.score, &its);

      if (bestIts)
	LogThisAt(3,"Deletion between " << pos-1 << " and " << pos << ": " << baseToChar(bestIts->base) << endl);

    } else if (isTMutStateIndex(mutState)) {

      const Pos dupIdx = tMutStateDupIdx (mutState);
      if (dupIdx < mdl - 1)
	updateBest (state, pos-1, tMutStateIndex(dupIdx+1), mutatorScores.sub[tanDupBase(ss,dupIdx+1)][seq[pos-1]], NULL);
      updateBest (state, pos, sMutStateIndex(), mutatorScores.tanDup + mutatorScores.len[dupIdx], NULL);

      if (bestMutState == sMutStateIndex()) {
	string dupstr;
	for (Pos dupIdx = tMutStateDupIdx(mutState); dupIdx >= 0; --dupIdx)
	  dupstr += baseToChar (tanDupBase(ss,dupIdx));
	LogThisAt(3,"Duplication at " << pos << ": " << dupstr << endl);
      }

    } else
      Abort ("Unknown traceback state");

    checkBest();
    if (bestInSym)
      trace.push_front (bestInSym);
  }

  return string (trace.begin(), trace.end());
}

vguard<FastSeq> decodeFastSeqs (const char* filename, const Machine& machine, const MutatorParams& mutatorParams) {
  const vguard<FastSeq> outseqs = readFastSeqs (filename);
  vguard<FastSeq> inseqs;
  InputModel inmod;
  for (auto& outseq: outseqs) {
    ViterbiMatrix vit (machine, inmod, mutatorParams, outseq);
    FastSeq inseq;
    inseq.name = outseq.name;
    inseq.seq = vit.traceback();
    inseqs.push_back (inseq);
  }
  return inseqs;
}
