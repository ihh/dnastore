#include <iomanip>
#include "builder.h"

vguard<int> TransBuilder::edgeFlagsToCountLookup ({ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 });

TransBuilder::TransBuilder (Pos len)
  : len (len),
    maxKmer (kmerMask (len)),
    maxTandemRepeatLen (len / 2),
    invertedRepeatLen (0),
    keepDegenerates (true),
    nControlWords (0),
    controlWordAtStart (false),
    controlWordAtEnd (false),
    startAndEndUseSameControlWord (false),
    buildDelayedMachine (false),
    kmerValid (maxKmer + 1)
{ }

void TransBuilder::findCandidates() {
  ProgressLog (plogReps, 1);
  plogReps.initProgress ("Filtering %d-mer repeats", len);
  kmers.clear();
  for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
    plogReps.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      
    if (!endsWithMotif(kmer,len,excludedMotif,"excluded motif")
	&& !endsWithMotif(kmer,len,excludedMotifRevComp,"revcomp of excluded motif")
	&& !hasExactTandemRepeat(kmer,len,maxTandemRepeatLen)
	&& !hasExactLocalInvertedRepeat(kmer,len,2,maxTandemRepeatLen)
	&& !hasExactNonlocalInvertedRepeat(kmer,len,invertedRepeatLen,2)) {
      LogThisAt(9,"Accepting " << kmerString(kmer,len) << endl);
      kmerValid[kmer] = true;
      kmers.push_back (kmer);
    }
  }
  const auto nKmersWithoutReps = kmers.size();
  LogThisAt(2,"Found " << nKmersWithoutReps << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutReps/(1.+(double)maxKmer) << "%)" << endl);
}

void TransBuilder::pruneUnreachable() {
  map<Kmer,Pos> dist;
  for (const auto& kl: sourceMotif)
    if (kl.len == len)
      doDFS (kl.kmer, dist);
  if (kmers.size() && !dist.size())
    doDFS (kmers.front(), dist);
  unsigned long long nDropped = 0;
  for (auto kmer: kmers)
    if (!dist.count(kmer)) {
      LogThisAt(6,"Dropping " << kmerString(kmer,len) << " as it was not seen in depth-first search" << endl);
      kmerValid[kmer] = false;
      ++nDropped;
    }
  if (nDropped) {
    LogThisAt(4,"Dropped " << nDropped << " " << len << "-mers that were unreachable in depth-first search" << endl);
    kmers.clear();
    for (const auto& kd: dist)
      kmers.push_back (kd.first);
    pruneDeadEnds();
  } else
    LogThisAt(5,"All " << kmers.size() << " " << len << "-mers were reached in depth-first search" << endl);
}

void TransBuilder::doDFS (Kmer kmer, map<Kmer,Pos>& distance) const {
  EdgeVector nbr;
  list<Kmer> kqueue;
  list<Pos> kdist;
  kqueue.push_back (kmer);
  kdist.push_back (0);
  while (!kqueue.empty()) {
    kmer = kqueue.back();
    const auto d = kdist.back();
    LogThisAt(9,"Depth-first search: visiting " << kmerString(kmer,len) << endl);
    kqueue.pop_back();
    kdist.pop_back();
    if (!distance.count(kmer)) {
      distance[kmer] = d;
      getOutgoing (kmer, nbr);
      for (auto n: nbr)
	if (kmerValid[n] && !distance.count(n)) {
	  kqueue.push_back (n);
	  kdist.push_back (d + 1);
	}
    }
  }
}

set<Kmer> TransBuilder::kmersEndingWith (KmerLen motif) const {
  set<Kmer> result;
  for (auto kmer: kmers)
    if (endsWithMotif(kmer,len,motif))
      result.insert(kmer);
  return result;
}

Pos TransBuilder::stepsToReach (KmerLen motif, int maxSteps) const {
  set<Kmer> nbr = kmersEndingWith(motif);
  EdgeVector in;
  for (int steps = 0; steps < maxSteps; ++steps) {
    if (nbr.size() == kmers.size())
      return steps;
    set<Kmer> prev;
    for (auto kmer: nbr)
      if (!(endsWithMotif(kmer,len,sourceMotif) || endsWithMotif(kmer,len,motif))
	  || steps == 0) {
	getIncoming (kmer, in);
	for (auto p: in)
	  if (kmerValid[p])
	    prev.insert (p);
      }
    nbr.swap (prev);
  }
  return -1;
}

map<Kmer,list<Kmer> > TransBuilder::pathsTo (Kmer dest, int steps) const {
  map<Kmer,list<Kmer> > pathFrom;
  pathFrom[dest] = list<Kmer>();
  EdgeVector in;
  for (int step = steps - 1; step >= 0; --step) {
    map<Kmer,list<Kmer> > longerPathFrom;
    for (auto interPath: pathFrom) {
      const Kmer inter = interPath.first;
      getIncoming (inter, in);
      for (auto src: in)
	if (kmerValid[src]
	    && (!(endsWithMotif(src,len,sourceMotif) || src == dest)
		|| step == 0))
	  (longerPathFrom[src] = interPath.second).push_front (inter);
    }
    pathFrom.swap (longerPathFrom);
  }
  return pathFrom;
}

void TransBuilder::pruneDeadEnds() {
  ProgressLog (plogPrune, 3);
  plogPrune.initProgress ("Pruning dead ends");
  for (auto kmer: kmers)
    pruneDeadEnds (kmer);
  const unsigned long long nKmers = kmers.size();
  unsigned long long nPruned = 0, nUnpruned = 0;
  list<Kmer> unprunedKmers;
  for (auto kmer: kmers) {
    ++nPruned;
    plogPrune.logProgress (nPruned / (double) nKmers, "sequence %llu/%llu", nPruned, nKmers);
    if (kmerValid[kmer]) {
      unprunedKmers.push_back (kmer);
      ++nUnpruned;
    }
  }
  LogThisAt(4,"Dead-end pruning removed " << (nPruned - nUnpruned) << " " << len << "-mers, leaving " << nUnpruned << endl);
  kmers.swap (unprunedKmers);
}

void TransBuilder::assertKmersCorrect() const {
  const set<Kmer> kmerSet (kmers.begin(), kmers.end());
  for (Kmer kmer = 0; kmer <= maxKmer; ++kmer)
    if (kmerValid[kmer])
      Assert (kmerSet.count(kmer), "Missing kmer %s from kmer list", kmerString(kmer,len).c_str());
  for (Kmer kmer: kmers)
    Assert (kmerValid[kmer], "Invalid kmer %s in kmer list", kmerString(kmer,len).c_str());
}

void TransBuilder::buildEdges() {
  LogThisAt(1,"Building edge graph for " << kmers.size() << " " << len << "-mers" << endl);
  EdgeVector out;
  for (auto kmer: kmers) {
    EdgeFlags outFlags = outgoingEdgeFlags(kmer,out);
    if (edgeFlagsToCount(outFlags) > 2 && !keepDegenerates) {
      if ((outFlags & PurineFlags) == PurineFlags)
	outFlags = dropWorseEdge (kmer, outFlags, out, AdenineBase, GuanineBase);
      if ((outFlags & PyrimidineFlags) == PyrimidineFlags)
	outFlags = dropWorseEdge (kmer, outFlags, out, CytosineBase, ThymineBase);
    }
    kmerOutFlags[kmer] = outFlags;
  }
  if (!keepDegenerates)
    LogThisAt(2,"Dropped " << droppedEdge.size() << " degenerate transitions" << endl);
  pruneDeadEnds();
}

void TransBuilder::indexStates() {
  assertKmersCorrect();
  
  nStates = 0;
  if (controlWordAtStart)
    nStates += buildDelayedMachine ? (len/2) : len;
  else
    nStates++;
  if (isStartControlIndex(0) && isEndControlIndex(0))
    ++nStates;
  for (auto kmer: controlWord)
    kmerState[kmer] = nStates++;
  for (auto kmer: kmers)
    if (!kmerState.count(kmer) && endsWithMotif(kmer,len,sourceMotif))
      kmerState[kmer] = nStates++;
  firstNonControlState = nStates;
  for (auto kmer: kmers)
    if (!kmerState.count(kmer))
      kmerState[kmer] = nStates++;
  for (auto kmer: kmers) {
    const auto nOut = countOutgoing(kmer);
    if (nOut > 2)
      kmerStateZero[kmer] = nStates++;
    if (nOut > 3)
      kmerStateOne[kmer] = nStates++;
  }
  for (size_t c = 0; c < nControlWords; ++c) {
    vguard<map<Kmer,State> > ckState (controlWordSteps[c]);
    for (Pos step = 0; step < controlWordSteps[c] - 1; ++step)
      for (auto kmer: controlWordIntermediates[c][step])
	ckState[step][kmer] = nStates++;
    controlKmerState.push_back (ckState);
  }

  if (buildDelayedMachine)
    nStates += len / 2;

  endState = nStates++;
}

void TransBuilder::prepare() {
  findCandidates();
  pruneDeadEnds();
  pruneUnreachable();
  getControlWords();
  buildEdges();
  indexStates();
}

Machine TransBuilder::makeMachine() {
  if (buildDelayedMachine) {
    Require (len % 2 == 0, "Delayed machine must have even number of bases per word");
    Require (controlWordAtStart && controlWordAtEnd && nControlWords > 0, "Delayed machine must generate control words at start & end of encoded sequence");
  }

  prepare();

  Machine machine;
  machine.state = vguard<MachineState> (nStates);

  if (controlWordAtStart) {
    const Pos p0 = buildDelayedMachine ? len/2 : 0;
    for (Pos p = p0; p < len; ++p) {
      const State s = p - p0;
      MachineState& ms = machine.state[s];
      ms.leftContext = string(len-p,'*') + kmerSubstring(startControlWord(),len-p+1,p);
      ms.name = (s == 0 ? "Start#" : "Load(Start)#") + to_string(s);
      ms.trans.push_back (MachineTransition (MachineNull, baseToChar(getBase(startControlWord(),len-p)), s+1));
    }
  } else {
    MachineState& ms = machine.state.front();
    ms.leftContext = string(len,'*');
    ms.name = "Start#1";
    ms.trans.push_back (MachineTransition (MachineNull, MachineNull, firstNonControlState));
  }
  
  EdgeVector out;
  vguard<char> outChar;
  vguard<State> outState;

  int nOut2 = 0, nOut3 = 0, nOut4 = 0;
  for (auto kmer: kmers) {
    const State s = kmerState.at(kmer);
    MachineState& ms = machine.state[s];
    ms.leftContext = kmerString(kmer,len);
    
    getOutgoing (kmer, out);
    const EdgeFlags outFlags = kmerOutFlags.at(kmer);
    outChar.clear();
    outState.clear();
    for (size_t n = 0; n < 4; ++n)
      if (outFlags & (1 << n)) {
	outChar.push_back (baseToChar(n));
	outState.push_back (kmerState.at(out[n]));
      }

    ms.name = "Code";
    if (endsWithMotif(kmer,len,sourceMotif)) {
      ms.name = "Source";
      for (size_t c = 0; c < controlWord.size(); ++c)
	if (kmer == controlWord[c]) {
	  if (isEndControlIndex(c))
	    ms.name = string("Control(End)");
	  else if (isStartControlIndex(c))
	    ms.name = string("Control(Start)");
	  else
	    ms.name = string("Control(") + controlChar(c) + ")";
	}
    }
    ms.name += "#" + to_string(s);

    if (outChar.size() == 1)
      ms.trans.push_back (MachineTransition (MachineNull, outChar[0], outState[0]));

    else if (outChar.size() == 2) {
      const int rotate2 = (++nOut2 % 2);
      const size_t i2 = rotate2, j2 = (rotate2 + 1) % 2;
      ms.trans.push_back (MachineTransition (MachineBit0, outChar[i2], outState[i2]));
      ms.trans.push_back (MachineTransition (MachineBit1, outChar[j2], outState[j2]));

      ms.trans.push_back (MachineTransition (MachineFlush, MachineNull, s));

      ms.trans.push_back (MachineTransition (MachineStrictBit0, outChar[i2], outState[i2]));
      ms.trans.push_back (MachineTransition (MachineStrictBit1, outChar[j2], outState[j2]));

    } else if (outChar.size() == 3) {
      const int rotate3 = (++nOut3 % 3);
      const size_t i3 = rotate3, j3 = (rotate3 + 1) % 3, k3 = (rotate3 + 2) % 3;
      const State s0 = kmerStateZero.at(kmer);
      ms.trans.push_back (MachineTransition (MachineBit0, MachineNull, s0));
      ms.trans.push_back (MachineTransition (MachineBit1, outChar[k3], outState[k3]));

      machine.state[s0].leftContext = kmerString(kmer,len);
      machine.state[s0].name = string("Split0#") + to_string(s0);
      machine.state[s0].trans.push_back (MachineTransition (MachineBit0, outChar[i3], outState[i3]));
      machine.state[s0].trans.push_back (MachineTransition (MachineBit1, outChar[j3], outState[j3]));

      ms.trans.push_back (MachineTransition (MachineFlush, MachineNull, s));
      machine.state[s0].trans.push_back (MachineTransition (MachineFlush, outChar[i3], outState[i3]));

      ms.trans.push_back (MachineTransition (MachineStrictTrit0, outChar[i3], outState[i3]));
      ms.trans.push_back (MachineTransition (MachineStrictTrit1, outChar[j3], outState[j3]));
      ms.trans.push_back (MachineTransition (MachineStrictTrit2, outChar[k3], outState[k3]));

    } else if (outChar.size() == 4) {
      const int rotate4 = (++nOut4 % 4);
      const size_t i4 = rotate4, j4 = (rotate4 + 1) % 4, k4 = (rotate4 + 2) % 4, l4 = (rotate4 + 3) % 4;
      const State s0 = kmerStateZero.at(kmer);
      const State s1 = kmerStateOne.at(kmer);
      ms.trans.push_back (MachineTransition (MachineBit0, MachineNull, s0));
      ms.trans.push_back (MachineTransition (MachineBit1, MachineNull, s1));

      machine.state[s0].leftContext = kmerString(kmer,len);
      machine.state[s0].name = string("Split0#") + to_string(s0);
      machine.state[s0].trans.push_back (MachineTransition (MachineBit0, outChar[i4], outState[i4]));
      machine.state[s0].trans.push_back (MachineTransition (MachineBit1, outChar[j4], outState[j4]));

      machine.state[s1].leftContext = kmerString(kmer,len);
      machine.state[s1].name = string("Split1#") + to_string(s1);
      machine.state[s1].trans.push_back (MachineTransition (MachineBit0, outChar[k4], outState[k4]));
      machine.state[s1].trans.push_back (MachineTransition (MachineBit1, outChar[l4], outState[l4]));

      ms.trans.push_back (MachineTransition (MachineFlush, MachineNull, s));
      machine.state[s0].trans.push_back (MachineTransition (MachineFlush, outChar[i4], outState[i4]));
      machine.state[s1].trans.push_back (MachineTransition (MachineFlush, outChar[l4], outState[l4]));

      ms.trans.push_back (MachineTransition (MachineStrictQuat0, outChar[i4], outState[i4]));
      ms.trans.push_back (MachineTransition (MachineStrictQuat1, outChar[j4], outState[j4]));
      ms.trans.push_back (MachineTransition (MachineStrictQuat2, outChar[k4], outState[k4]));
      ms.trans.push_back (MachineTransition (MachineStrictQuat3, outChar[l4], outState[l4]));
    }
    
    if (outChar.size() > 1) {
      for (size_t c = 0; c < controlWord.size(); ++c) {
	if (isSourceControlIndex(c))
	  continue;
	ms.trans.push_back (controlTrans (s, controlWordPath[c].at(kmer).front(), c, 0));
      }
      if (!controlWordAtEnd)
	ms.trans.push_back (MachineTransition (MachineEOF, 0, endState));
    }
    if (controlWordAtEnd && kmer == endControlWord()) {

      if (isStartControlIndex(0) && isEndControlIndex(0)) {
	swap (ms.trans, machine.state[s-1].trans);  // give the Start copy all the outgoing transitions
	machine.state[s-1].leftContext = ms.leftContext;
	machine.state[s-1].name = string("Control(Start)") + "#" + to_string(s-1);
      }

      if (buildDelayedMachine)
	ms.trans.push_back (MachineTransition (0, '*', endState - len/2));
      else
	ms.trans.push_back (MachineTransition (0, 0, endState));
    }
  }

  for (size_t c = 0; c < controlWord.size(); ++c)
    for (int step = 0; step < controlWordSteps[c] - 1; ++step) {
      const auto& ckState = controlKmerState[c][step];
      for (const auto& ks: ckState) {
	const Kmer srcKmer = ks.first;
	const State srcState = ks.second;
	const Kmer destKmer = nextIntermediateKmer (srcKmer, c, step + 1);
	machine.state[srcState].leftContext = kmerString(srcKmer,len);
	machine.state[srcState].name = (isEndControlIndex(c) ? string("Bridge(End)") : (string("Bridge(") + controlChar(c) + ")")) + "#" + to_string(srcState);
	machine.state[srcState].trans.push_back (controlTrans (srcState, destKmer, c, step + 1));
      }
    }

  machine.state[endState].name = "End#" + to_string(endState);
  machine.state[endState].leftContext = buildDelayedMachine
    ? (kmerSubstring(endControlWord(),1,len/2) + string(len/2,'*'))
    : (controlWordAtEnd ? kmerString(endControlWord(),len) : string(len,'*'));
  
  if (buildDelayedMachine) {
    for (Pos pos = 1; pos <= len/2; ++pos) {
      const State s = endState - 1 - len/2 + pos;
      MachineState& ms = machine.state[s];
      ms.name = string("Unload(End)#") + to_string(s);
      ms.leftContext = kmerSubstring(endControlWord(),1,len-pos) + string(pos,'*');
      if (pos < len/2)
	ms.trans.push_back (MachineTransition (0, '*', s + 1));
      else
	ms.trans.push_back (MachineTransition (0, 0, endState));
    }

    for (auto& ms: machine.state) {
      ms.rightContext = string (ms.leftContext.begin() + len/2, ms.leftContext.end());
      ms.leftContext.erase (ms.leftContext.begin() + len/2, ms.leftContext.end());
    }

    for (auto& ms: machine.state)
      for (auto& t: ms.trans)
	if (t.out)
	  t.out = machine.state[t.dest].leftContext.back();
  }
  
  return machine;
}

bool TransBuilder::isSourceControlIndex (size_t c) const {
  return isStartControlIndex(c) && !isEndControlIndex(c);
}

bool TransBuilder::isStartControlIndex (size_t c) const {
  return controlWordAtStart && c == 0;
}

bool TransBuilder::isEndControlIndex (size_t c) const {
  return controlWordAtEnd && c == (controlWordAtStart && !startAndEndUseSameControlWord ? 1 : 0);
}

Kmer TransBuilder::startControlWord() const {
  Assert (controlWordAtStart, "There is no start control word");
  return controlWord.front();
}

Kmer TransBuilder::endControlWord() const {
  Assert (controlWordAtEnd, "There is no end control word");
  return controlWord[controlWordAtStart && !startAndEndUseSameControlWord ? 1 : 0];
}

char TransBuilder::controlChar (size_t nControlWord) const {
  const size_t nReserved = (controlWordAtStart ? 1 : 0) + (controlWordAtEnd ? 1 : 0);
  return Machine::controlChar (nControlWord - nReserved);
}

MachineTransition TransBuilder::controlTrans (State srcState, Kmer destKmer, size_t nControlWord, size_t step) const {
  const State destState =
    (step == controlWordSteps[nControlWord] - 1 && destKmer == controlWord[nControlWord])
    ? kmerState.at(destKmer)
    : controlKmerState[nControlWord][step].at(destKmer);
  return MachineTransition (step == 0
			    ? (isEndControlIndex(nControlWord) ? MachineEOF : controlChar(nControlWord))
			    : MachineNull, baseToChar(getBase(destKmer,1)), destState);
}

Kmer TransBuilder::nextIntermediateKmer (Kmer srcKmer, size_t nControlWord, size_t step) const {
  EdgeVector out;
  getOutgoing (srcKmer, out);
  for (Kmer destKmer: out)
    if ((step == controlWordSteps[nControlWord] - 1 && destKmer == controlWord[nControlWord])
	|| (step < controlWordSteps[nControlWord] - 1 && controlWordIntermediates[nControlWord][step].count(destKmer)))
      return destKmer;
  Abort("Can't find intermediate kmer following %s at step %d to control word #%d (%s)", kmerString(srcKmer,len).c_str(), step, nControlWord, kmerString(controlWord[nControlWord],len).c_str());
  return 0;
}

bool TransBuilder::getNextControlWord() {
  if (controlWord.size() == nControlWords)
    return true;
  const size_t cCurrent = controlWord.size();
  const bool currentIsSource = isSourceControlIndex(cCurrent);
  LogThisAt(3,"Looking for control word #" << (cCurrent + 1)
	    << (controlWord.empty() ? string() : (string(" (previous: ") + to_string_join(controlWordString) + ")"))
	    << endl);
  vguard<Kmer> cand (kmers.begin(), kmers.end());
  vguard<size_t> dist (cand.size(), len);
  for (size_t k = 0; k < cand.size(); ++k)
    if (kmerValid[cand[k]])
      for (auto cw: controlWord)
	dist[k] = min (dist[k], min (kmerHammingDistance (cand[k], cw, len),
				     kmerHammingDistance (cand[k], kmerRevComp(cw,len), len)));
  auto indexByDistance = orderedIndices (dist);
  while (!indexByDistance.empty()) {
    const size_t bestIdx = indexByDistance.back();
    indexByDistance.pop_back();
    if (dist[bestIdx] == 0)
      continue;
    const Kmer best = cand[bestIdx];
    const KmerLen bestMotif (best, len);
    const Pos steps = stepsToReach (bestMotif);
    if (!currentIsSource && steps < 0) {
      LogThisAt(5,"Rejecting " << kmerString(bestMotif)
		<< " for control word #" << (cCurrent + 1)
		<< " as it is not reachable" << endl);
      continue;
    }

    const Kmer bestRevComp = kmerRevComp (best, len);
    if (bestRevComp == best) {
      LogThisAt(5,"Rejecting " << kmerString(bestMotif)
		<< " for control word #" << (cCurrent + 1)
		<< " as it is palindromic" << endl);
      continue;
    }

    LogThisAt(3,"Trying control word " << kmerString(bestMotif)
	      << (currentIsSource ? string() : (string(" which is reachable in ") + to_string(steps) + " steps"))
	      << (controlWord.empty()
		  ? string()
		    : (string(" and has ")
		       + to_string(dist[bestIdx])
		       + "+ differences from (" + to_string_join(controlWordString) + ")"))
	      << endl);

    list<Kmer> savedKmers;
    for (auto kmer: kmers)
      if (kmerValid[kmer])
	savedKmers.push_back (kmer);
    
    sourceMotif.insert (bestMotif);
    kmerValid[bestRevComp] = false;

    pruneDeadEnds();
    pruneUnreachable();

    bool broken = false;
    if (!currentIsSource && stepsToReach (bestMotif) < 0) {
      LogThisAt(4,"Oops - " << kmerString(bestMotif) << " is unreachable when reverse-complement " << kmerString(bestRevComp,len) << " is excluded" << endl);
      broken = true;
    }
    
    for (size_t c = 0; !broken && c < cCurrent; ++c)
      if (!isSourceControlIndex(c)) {
	const Kmer prevControl = controlWord[c];
	const KmerLen prevMotif (prevControl, len);
	const Pos prevSteps = stepsToReach (prevMotif);
	if (prevSteps < 0) {
	  LogThisAt(4,"Oops - setting " << kmerString(bestMotif) << " as a control word breaks paths to previous control word " << kmerString(prevMotif) << endl);
	  broken = true;
	}
      }

    if (!broken) {
      controlWord.push_back (best);
      controlWordString.push_back (kmerString(best,len));

      if (getNextControlWord())
	return true;
    
      controlWord.pop_back();
      controlWordString.pop_back();
    }

    // flag this word as unusable and restore previous state
    dist[bestIdx] = 0;
    sourceMotif.erase (bestMotif);
    for (auto kmer: savedKmers)
      kmerValid[kmer] = true;

    LogThisAt(3,"Trying next option for control word #" << (cCurrent + 1) << endl);
  }
  return false;
}

void TransBuilder::getControlWords() {
  if (nControlWords > 0)
    LogThisAt(1,"Attempting to allocate " << plural(nControlWords,"control word") << endl);
  
  if (nControlWords == 0 && controlWordAtStart) {
    Warn ("No control words allocated, disabling control word at start");
    controlWordAtStart = false;
  }

  if (nControlWords == 0 && controlWordAtEnd) {
    Warn ("No control words allocated, disabling control word at end");
    controlWordAtEnd = false;
  }

  if (nControlWords == 1 && controlWordAtStart && controlWordAtEnd && !startAndEndUseSameControlWord) {
    Warn ("Only 1 control word allocated, so start and end will use same control word");
    startAndEndUseSameControlWord = true;
  }
  
  Require (getNextControlWord(), "Ran out of control words");

  if (controlWordAtEnd && (!startAndEndUseSameControlWord || !controlWordAtStart)) {
    const Kmer e = endControlWord();
    EdgeVector out;
    getOutgoing (e, out);
    for (auto n: out)
      if (kmerValid[n])
	droppedEdge.insert (pair<Kmer,Kmer> (e, n));
  }
  
  pruneDeadEnds();
  pruneUnreachable();

  size_t totalInter = 0;
  for (size_t c = 0; c < controlWord.size(); ++c) {
    const Kmer controlKmer = controlWord[c];
    if (isSourceControlIndex(c)) {
      controlWordSteps.push_back (0);
      controlWordPath.push_back (map<Kmer,list<Kmer> >());
      controlWordIntermediates.push_back (vguard<set<Kmer> >());
    } else {
      const Pos controlSteps = stepsToReach (KmerLen (controlWord[c], len));
      Assert (controlSteps >= 0, "Control word #%d unreachable", c+1);
      controlWordSteps.push_back (controlSteps);
      controlWordPath.push_back (pathsTo (controlKmer, controlSteps));
      vguard<set<Kmer> > intermediates (controlSteps);
      for (auto kmer: kmers)
	if (!controlWordAtEnd || kmer != endControlWord() || (startAndEndUseSameControlWord && controlWordAtStart)) {
	  Pos step = 0;
	  for (auto inter: controlWordPath[c].at(kmer)) {
	    LogThisAt(9,"Adding " << kmerString(inter,len) << " at step " << step << " from " << kmerString(kmer,len) << " to control word #" << c << " (" << kmerString(controlWord[c],len) << ")" << endl);
	    intermediates[step++].insert (inter);
	  }
	}
      intermediates.pop_back();
      size_t nInter = 0;
      for (const auto& ks: intermediates)
	nInter += ks.size();
      controlWordIntermediates.push_back (intermediates);
      LogThisAt(3,"Control word " << kmerString(controlKmer,len) << " needs " << nInter << " intermediate states" << endl);
      totalInter += nInter;
    }
  }
  LogThisAt(2,"Control words (" << join(controlWordString) << ") require " << totalInter << " bridge states" << endl);
}
