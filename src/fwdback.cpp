#include "fwdback.h"

MutatorMatrix::MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : mutatorParams (mutatorParams),
    mutatorScores (mutatorParams),
    stock (stock),
    cell (nCells (mutatorParams, stock)),
    maxDupLen (mutatorParams.maxDupLen()),
    alignCols (stock.columns()),
    align (stock.gapped)
{
  Assert (stock.rows() == 2, "Training mutator model requires a 2-row alignment; this alignment has %d rows", stock.rows());
  const string dna (dnaAlphabet);
  inSeq = align.ungapped[0].tokens(dna);
  outSeq = align.ungapped[1].tokens(dna);
  const AlignRowPath& inPath = align.path.at(0);
  const AlignRowPath& outPath = align.path.at(1);
  SeqIdx inSeqIdx = 0, outSeqIdx = 0;
  for (AlignColIndex col = 0; col < alignCols; ++col) {
    if (inPath[col])
      ++inSeqIdx;
    if (outPath[col])
      ++outSeqIdx;
    col2InSeqIdx.push_back (inSeqIdx);
    col2OutSeqIdx.push_back (outSeqIdx);
  }
}

FwdBackMatrix::FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : fwd (mutatorParams, stock),
    back (mutatorParams, stock)
{
  // WRITE ME
}

MutatorCounts FwdBackMatrix::counts() const {
  MutatorCounts mc (fwd.mutatorParams);
  // WRITE ME
  return mc;
}

    
