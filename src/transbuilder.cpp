#include <iomanip>
#include "transbuilder.h"

vguard<int> TransBuilder::edgeFlagsToCountLookup ({ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 });

TransBuilder::TransBuilder (Pos len)
  : len (len),
    maxKmer (kmerMask (len)),
    maxTandemRepeatLen (len / 2),
    invertedRepeatLen (0),
    controlWords (0),
    kmerValid (maxKmer + 1)
{ }

void TransBuilder::removeRepeats() {
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
  map<Kmer,int> dist;
  for (const auto& kl: sourceMotif)
    if (kl.len == len)
      doDFS (kl.kmer, dist, true);
  if (kmers.size() && !dist.size())
    doDFS (kmers.front(), dist, true);
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

void TransBuilder::doDFS (Kmer kmer, map<Kmer,int>& distance, bool forwards) const {
  EdgeVector nbr;
  list<Kmer> kqueue;
  list<int> kdist;
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
      if (forwards)
	getOutgoing (kmer, nbr);
      else
	getIncoming (kmer, nbr);
      for (auto n: nbr)
	if (kmerValid[n] && !distance.count(n)) {
	  kqueue.push_back (n);
	  kdist.push_back (d + 1);
	}
    }
  }
}

set<Kmer> TransBuilder::neighbors (const set<Kmer>& start, int steps, bool forwards) const {
  if (steps > 1)
    return neighbors (neighbors(start,steps-1,forwards), 1, forwards);
  set<Kmer> nbr;
  EdgeVector next;
  for (auto kmer: start) {
    if (forwards)
      getOutgoing (kmer, next);
    else
      getIncoming (kmer, next);
    for (auto n: next)
      if (kmerValid[n])
	nbr.insert (n);
  }
  return nbr;
}

set<Kmer> TransBuilder::neighbors (Kmer start, int steps, bool forwards) const {
  set<Kmer> startSet;
  startSet.insert (start);
  return neighbors (startSet, steps, forwards);
}

int TransBuilder::stepsToReach (KmerLen motif, int maxSteps) const {
  set<KmerLen> motifs;
  motifs.insert (motif);
  set<Kmer> nbr;
  for (auto kmer: kmers)
    if (endsWithMotif(kmer,len,motifs))
      nbr.insert(kmer);
  for (int steps = 0; steps < maxSteps; ++steps) {
    if (nbr.size() == kmers.size())
      return steps;
    nbr = neighbors (nbr, 1, false);
  }
  return -1;
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
    EdgeFlags outFlags = outgoingEdgeFlags(kmer,out);
    // WRITE ME: if we have a transition degeneracy, skip one of them
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
    // number of steps with which this state can be reached
    //    outs << " (" << stepsToReach(KmerLen(kmer,len),len*2) << ")";
    outs << endl;
  }
}


