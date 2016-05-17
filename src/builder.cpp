#include <iomanip>
#include "builder.h"

vguard<int> TransBuilder::edgeFlagsToCountLookup ({ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 });

TransBuilder::TransBuilder (Pos len)
  : len (len),
    maxKmer (kmerMask (len)),
    maxTandemRepeatLen (len / 2),
    invertedRepeatLen (0),
    keepDegenerates (false),
    nControlWords (0),
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
	&& !hasExactLocalInvertedRepeat(kmer,len,3,maxTandemRepeatLen)
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
    LogThisAt(2,"Dropped " << nDropped << " " << len << "-mers that were unreachable in depth-first search" << endl);
    kmers.clear();
    for (const auto& kd: dist)
      kmers.push_back (kd.first);
    pruneDeadEnds();
  } else
    LogThisAt(2,"All " << kmers.size() << " " << len << "-mers were reached in depth-first search" << endl);
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
  const size_t nKmersWithoutMotif = kmers.size() - nbr.size();
  EdgeVector in;
  for (int steps = 0; steps < maxSteps; ++steps) {
    if (nbr.size() == nKmersWithoutMotif)
      return steps;
    set<Kmer> prev;
    for (auto kmer: nbr) {
      getIncoming (kmer, in);
      for (auto p: in)
	if (kmerValid[p] && (!endsWithMotif(p,len,motif) || steps == maxSteps-2))
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
  for (int step = 0; step < steps; ++step) {
    map<Kmer,list<Kmer> > longerPathFrom;
    for (auto interPath: pathFrom) {
      const Kmer inter = interPath.first;
      getIncoming (inter, in);
      for (auto src: in)
	if (kmerValid[src] && (src != dest || step == steps - 1))
	  (longerPathFrom[src] = interPath.second).push_front (inter);
    }
    pathFrom.swap (longerPathFrom);
  }
  return pathFrom;
}

void TransBuilder::pruneDeadEnds() {
  ProgressLog (plogPrune, 1);
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
  LogThisAt(2,"Dead-end pruning removed " << (nPruned - nUnpruned) << " " << len << "-mers, leaving " << nUnpruned << endl);
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
  EdgeVector out;
  for (auto kmer: kmers) {
    EdgeFlags outFlags = outgoingEdgeFlags(kmer,out);
    if (!keepDegenerates) {
      if ((outFlags & 3) == 3)  // A,G
	outFlags = dropWorseEdge (kmer, outFlags, out, 0, 1);
      if ((outFlags & 12) == 12)  // C,T
	outFlags = dropWorseEdge (kmer, outFlags, out, 2, 3);
    }
    kmerOutFlags[kmer] = outFlags;
  }
  if (!keepDegenerates)
    LogThisAt(2,"Dropped " << droppedEdge.size() << " degenerate edges" << endl);
  pruneDeadEnds();
}

void TransBuilder::indexStates() {
  nStates = 0;
  for (auto kmer: kmers)
    kmerState[kmer] = nStates++;
  for (auto kmer: kmers) {
    const auto nOut = countOutgoing(kmer);
    if (nOut > 2)
      kmerStateZero[kmer] = nStates++;
    if (nOut > 3)
      kmerStateOne[kmer] = nStates++;
  }
  nCodingStates = nStates;
  for (size_t c = 0; c < nControlWords; ++c) {
    vguard<map<Kmer,State> > ckState (controlWordSteps[c]);
    for (Pos step = 0; step < controlWordSteps[c] - 1; ++step)
      for (auto kmer: controlWordIntermediates[c][step])
	ckState[step][kmer] = nStates++;
    controlKmerState.push_back (ckState);
  }
}

Machine TransBuilder::makeMachine() {
  Machine machine (len);
  machine.state = vguard<MachineState> (nStates);

  EdgeVector out;
  vguard<char> outChar;
  vguard<State> outState;
  
  for (auto kmer: kmers) {
    const State s = kmerState.at(kmer);
    MachineState& ms = machine.state[s];
    ms.context = kmer;
    ms.type = 0;
    
    getOutgoing (kmer, out);
    const EdgeFlags outFlags = kmerOutFlags.at(kmer);
    outChar.clear();
    outState.clear();
    for (size_t n = 0; n < 4; ++n)
      if (outFlags & (1 << n)) {
	outChar.push_back (baseToChar(n));
	outState.push_back (kmerState.at(out[n]));
      }

    if (outChar.size() == 1)
      ms.trans.push_back (MachineTransition ('\0', outChar[0], outState[0]));
    else if (outChar.size() == 2) {
      ms.trans.push_back (MachineTransition ('0', outChar[0], outState[0]));
      ms.trans.push_back (MachineTransition ('1', outChar[1], outState[1]));
    } else if (outChar.size() == 3) {
      const State s0 = kmerStateZero.at(kmer);
      ms.trans.push_back (MachineTransition ('0', '\0', s0));
      ms.trans.push_back (MachineTransition ('1', outChar[2], outState[2]));
      machine.state[s0].context = kmer;
      machine.state[s0].type = 0;
      machine.state[s0].trans.push_back (MachineTransition ('0', outChar[0], outState[0]));
      machine.state[s0].trans.push_back (MachineTransition ('1', outChar[1], outState[1]));
    } else if (outChar.size() == 4) {
      const State s0 = kmerStateZero.at(kmer);
      const State s1 = kmerStateOne.at(kmer);
      ms.trans.push_back (MachineTransition ('0', '\0', s0));
      ms.trans.push_back (MachineTransition ('1', '\0', s1));
      machine.state[s0].context = kmer;
      machine.state[s0].type = 0;
      machine.state[s0].trans.push_back (MachineTransition ('0', outChar[0], outState[0]));
      machine.state[s0].trans.push_back (MachineTransition ('1', outChar[1], outState[1]));
      machine.state[s1].context = kmer;
      machine.state[s1].type = 0;
      machine.state[s1].trans.push_back (MachineTransition ('0', outChar[2], outState[2]));
      machine.state[s1].trans.push_back (MachineTransition ('1', outChar[3], outState[3]));
    }
    
    if (outChar.size() > 1)
      for (size_t c = 0; c < controlWord.size(); ++c)
	ms.trans.push_back (controlTrans (s, controlWordPath[c].at(kmer).front(), c, 0));
  }

  for (size_t c = 0; c < controlWord.size(); ++c) {
    machine.control[controlChar(c)] = controlWord[c];
    machine.state[kmerState.at(controlWord[c])].type = controlChar(c);
    for (size_t step = 0; step < controlWordSteps[c] - 1; ++step) {
      const auto& ckState = controlKmerState[c][step];
      for (const auto& ks: ckState) {
	const Kmer srcKmer = ks.first;
	const State srcState = ks.second;
	const Kmer destKmer = nextIntermediateKmer (srcKmer, c, step + 1);
	machine.state[srcState].context = srcKmer;
	machine.state[srcState].type = controlChar(c);
	machine.state[srcState].trans.push_back (controlTrans (srcState, destKmer, c, step + 1));
      }
    }
  }

  return machine;
}

char TransBuilder::controlChar (size_t nControlWord) const {
  return 'a' + nControlWord;
}

MachineTransition TransBuilder::controlTrans (State srcState, Kmer destKmer, size_t nControlWord, size_t step) const {
  const State destState =
    (step == controlWordSteps[nControlWord] - 1 && destKmer == controlWord[nControlWord])
    ? kmerState.at(destKmer)
    : controlKmerState[nControlWord][step].at(destKmer);
  return MachineTransition (step == 0 ? controlChar(nControlWord) : '\0', baseToChar(getBase(destKmer,1)), destState);
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

void TransBuilder::getControlWords() {
  if (nControlWords == 0)
    return;
  vguard<Kmer> cand;
  vguard<Pos> steps;
  vguard<size_t> dist;
  ProgressLog (plogControl, 1);
  plogControl.initProgress ("Looking for control words");
  const size_t nKmers = kmers.size();
  size_t n = 0;
  for (auto kmer: kmers) {
    plogControl.logProgress (n / (double) nKmers, "sequence %llu/%llu", n, nKmers);
    ++n;
    const int s = stepsToReach(KmerLen(kmer,len));
    if (s >= 0) {
      cand.push_back (kmer);
      steps.push_back (s);
      dist.push_back (len);
    }
  }
  LogThisAt(4,"Found " << cand.size() << " potential control words" << endl);
  Require (nControlWords*2 < cand.size(), "Not enough control words");
  while (controlWord.size() < nControlWords) {
    size_t bestIdx = cand.size();
    for (size_t idx = 1; idx < cand.size(); ++idx)
      if (kmerValid[cand[idx]]
	  && (bestIdx == cand.size()
	      || dist[idx] > dist[bestIdx]
	      || (dist[idx] == dist[bestIdx] && steps[idx] < steps[bestIdx])))
	bestIdx = idx;
    Require (bestIdx < cand.size(), "Not enough control words");
    const Kmer best = cand[bestIdx];
    const Kmer bestRevComp = kmerRevComp (best, len);
    LogThisAt(3,"Selected control word " << kmerString(best,len)
	      << " which is reachable from any other words in " << steps[bestIdx] << " steps"
	      << (controlWord.empty()
		  ? string()
		    : (string(" and has at least ")
		       + to_string(dist[bestIdx])
		       + " bases different from all other control words"))
	      << endl);
    controlWord.push_back (best);
    controlWordSteps.push_back (steps[bestIdx]);
    controlWordPath.push_back (pathsTo (best, steps[bestIdx]));
    for (size_t k = 0; k < cand.size(); ++k)
      if (kmerValid[cand[k]])
	dist[k] = min (dist[k], min (kmerHammingDistance (cand[k], best, len),
				     kmerHammingDistance (cand[k], bestRevComp, len)));

    sourceMotif.insert (KmerLen (best, len));
    kmerValid[bestRevComp] = false;
  
    pruneDeadEnds();
    pruneUnreachable();
  }

  assertKmersCorrect();
  
  for (size_t c = 0; c < controlWord.size(); ++c) {
    const Kmer controlKmer = controlWord[c];
    vguard<set<Kmer> > intermediates (controlWordSteps[c]);
    for (auto kmer: kmers) {
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
  }
}
