#include <iomanip>
#include "transbuilder.h"

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
  set<KmerLen> motifs;
  motifs.insert (motif);
  for (auto kmer: kmers)
    if (endsWithMotif(kmer,len,motifs))
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
    for (auto kmer: nbr) {
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
  for (int step = 0; step < steps; ++step) {
    map<Kmer,list<Kmer> > longerPathFrom;
    for (auto destPath: pathFrom) {
      getIncoming (destPath.first, in);
      for (auto src: in)
	if (kmerValid[src])
	  if (!longerPathFrom.count (src))
	    (longerPathFrom[src] = destPath.second).push_front (destPath.first);
    }
    pathFrom.swap (longerPathFrom);
  }
  Assert (pathFrom.size() == kmers.size(), "Incomplete path map: got %llu kmers, expected %llu", pathFrom.size(), kmers.size());
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
  State n = 0;
  for (auto kmer: kmers)
    kmerState[kmer] = ++n;
}

void TransBuilder::output (ostream& outs) {
  EdgeVector out;
  vguard<char> outChar;
  vguard<State> outState;
  for (auto kmer: kmers) {
    getOutgoing (kmer, out);
    const EdgeFlags outFlags = kmerOutFlags.at(kmer);
    outChar.clear();
    outState.clear();
    for (size_t n = 0; n < 4; ++n)
      if (outFlags & (1 << n)) {
	outChar.push_back (baseToChar(n));
	outState.push_back (kmerState.at(out[n]));
      }
    outs << "#" << kmerState.at(kmer) << " " << kmerString(kmer,len);
    if (outChar.size() == 1)
      outs << " /" << outChar[0] << "->#" << outState[0];
    else if (outChar.size() == 2)
      outs << " 0/" << outChar[0] << "->#" << outState[0]
	   << " 1/" << outChar[1] << "->#" << outState[1];
    else if (outChar.size() == 3)
      outs << " 00/" << outChar[0] << "->#" << outState[0]
	   << " 01/" << outChar[1] << "->#" << outState[1]
	   << " 1/" << outChar[2] << "->#" << outState[2];
    else if (outChar.size() == 4)  // kind of silly/obvious, but leave it in for completeness
      outs << " 00/" << outChar[0] << "->#" << outState[0]
	   << " 01/" << outChar[1] << "->#" << outState[1]
	   << " 10/" << outChar[2] << "->#" << outState[2]
	   << " 11/" << outChar[3] << "->#" << outState[3];
    if (outChar.size() > 1)
      for (size_t c = 0; c < controlWord.size(); ++c) {
	const list<Kmer>& p = controlWordPath[c].at (kmer);
	outs << ' ' << (char) ('A' + c) << "/";
	for (auto step: p)
	  outs << baseToChar(getBase(step,1));
	outs << "->#CONTROL" << (c+1);
      }
    outs << endl;
  }
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
  if (nControlWords*2 < cand.size())
    while (controlWord.size() < nControlWords) {
      size_t bestIdx = 0;
      for (size_t idx = 1; idx < cand.size(); ++idx)
	if (dist[idx] > dist[bestIdx]
	    || (dist[idx] == dist[bestIdx] && steps[idx] < steps[bestIdx]))
	  bestIdx = idx;
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
	dist[k] = min (dist[k], min (kmerHammingDistance (cand[k], best, len),
				     kmerHammingDistance (cand[k], bestRevComp, len)));
    }

  for (auto cw: controlWord) {
    sourceMotif.insert (KmerLen (cw, len));
    kmerValid[kmerRevComp (cw, len)] = false;
  }
  
  pruneDeadEnds();
  pruneUnreachable();

  assertKmersCorrect();
  
  for (size_t c = 0; c < controlWord.size(); ++c) {
    const Kmer controlKmer = controlWord[c];
    set<Kmer> intermediates;
    for (auto kmer: kmers)
      for (auto inter: controlWordPath[c].at(kmer))
	if (inter != controlKmer)
	  intermediates.insert (inter);
    controlWordIntermediates.push_back (intermediates);
    LogThisAt(3,"Control word " << kmerString(controlKmer,len) << " needs " << intermediates.size() << " intermediate states" << endl);
  }
}
