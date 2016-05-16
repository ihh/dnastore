#include <iomanip>
#include "transbuilder.h"

vguard<int> TransBuilder::edgeFlagsToCountLookup ({ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 });

TransBuilder::TransBuilder (Pos len)
  : len (len),
    maxKmer (kmerMask (len)),
    kmerValid (maxKmer + 1),
    maxTandemRepeatLen (len / 2),
    invertedRepeatLen (0),
    keepDegenerates (false)
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

void TransBuilder::removeDegenerates() {
  if (keepDegenerates)
    return;
  list<Kmer> kmersWithoutTransitions;
  for (Kmer kmer: kmers) {
    const Kmer transKmer = makeTransition(kmer,1);
    if (kmerValid[kmer] && kmerValid[transKmer] && kmerEqualOrBetter(transKmer,kmer,len)) {
      kmerValid[kmer] = false;
      LogThisAt(4,"Eliminating "
		<< kmerString(kmer,len) << " (" << setprecision(2) << 100*gcContent(kmer,len) << "% GC, " << setw(3) << kmerEntropy(kmer,len) << " bits)"
		<< " since we also have "
		<< kmerString(transKmer,len) << " (" << setprecision(2) << 100*gcContent(transKmer,len) << "% GC, " << setw(3) << kmerEntropy(transKmer,len) << " bits)"
		<< endl);
    } else
      kmersWithoutTransitions.push_back (kmer);
    }
  const unsigned long long nKmersWithoutTransitions = kmersWithoutTransitions.size();
  LogThisAt(2,"Removed " << (kmers.size() - nKmersWithoutTransitions) << " transition redundancies leaving " << nKmersWithoutTransitions << " candidate " << len << "-mers without repeats (" << setprecision(2) << 100*(double)nKmersWithoutTransitions/(1.+(double)maxKmer) << "%)" << endl);
  kmers.swap (kmersWithoutTransitions);
}

void TransBuilder::pruneDeadEnds() {
  const Kmer maxKmer = kmerMask(len);
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
  LogThisAt(2,"Dead-end pruning left " << nUnpruned << " " << len << "-mers (" << setprecision(2) << 100*(double)nUnpruned/(1.+(double)maxKmer) << "%)" << endl);
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
    outs << endl;
  }
}


