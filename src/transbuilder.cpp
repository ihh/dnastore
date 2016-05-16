#include <iomanip>
#include "transbuilder.h"

vguard<int> TransBuilder::edgeFlagsToCountLookup ({ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 });

TransBuilder::TransBuilder (Pos len)
  : len (len),
    maxKmer (kmerMask (len)),
    kmerValid (maxKmer + 1),
    maxTandemRepeatLen (len / 2),
    invertedRepeatLen (0)
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

void TransBuilder::keepDFS() {
  set<Kmer> seen;
  for (const auto& kl: sourceMotif)
    if (kl.len == len)
      doDFS (kl.kmer, seen);
  if (kmers.size() && !seen.size())
    doDFS (kmers.front(), seen);
  unsigned long long nDropped = 0;
  for (auto kmer: kmers)
    if (!seen.count(kmer)) {
      LogThisAt(6,"Dropping " << kmerString(kmer,len) << " as it was not seen in depth-first search" << endl);
      kmerValid[kmer] = false;
      ++nDropped;
    }
  if (nDropped) {
    LogThisAt(2,"Dropped " << nDropped << " " << len << "-mers that were unreachable in depth-first search" << endl);
    list<Kmer> seenKmers (seen.begin(), seen.end());
    kmers.swap (seenKmers);
    pruneDeadEnds();
  } else
    LogThisAt(2,"All " << kmers.size() << " " << len << "-mers were reached in depth-first search" << endl);
}

void TransBuilder::doDFS (Kmer kmer, set<Kmer>& seen) const {
  vguard<Kmer> out (4);
  list<Kmer> kqueue;
  kqueue.push_back (kmer);
  while (!kqueue.empty()) {
    kmer = kqueue.back();
    LogThisAt(9,"Depth-first search: visiting " << kmerString(kmer,len) << endl);
    kqueue.pop_back();
    if (!seen.count(kmer)) {
      seen.insert (kmer);
      getOutgoing (kmer, out);
      for (auto dest: out)
	if (kmerValid[dest] && !seen.count(dest))
	  kqueue.push_back (dest);
    }
  }
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
  vguard<Kmer> out (4);
  vguard<char> outChar;
  vguard<State> outState;
  for (auto kmer: kmers) {
    const EdgeFlags outFlags = outgoingEdgeFlags(kmer,out);
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
    outs << endl;
  }
}


