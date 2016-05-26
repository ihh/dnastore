#include "fwdback.h"
#include "logsumexp.h"

#define BaumWelchMinFracInc .001
#define BaumWelchMaxIter 100

MutatorMatrix::MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : mutatorParams (mutatorParams),
    mutatorScores (mutatorParams),
    stock (stock),
    cell (nCells (mutatorParams, stock), -numeric_limits<double>::infinity()),
    maxDupLen (mutatorParams.maxDupLen()),
    alignCols (stock.columns()),
    align (stock.gapped),
    inPath (align.path.at(0)),
    outPath (align.path.at(1))
{
  Assert (stock.rows() == 2, "Training mutator model requires a 2-row alignment; this alignment has %d rows", stock.rows());
  const string dna (dnaAlphabet);
  inSeq = align.ungapped[0].tokens(dna);
  outSeq = align.ungapped[1].tokens(dna);
  SeqIdx inSeqIdx = 0, outSeqIdx = 0;
  for (AlignColIndex col = 0; col < alignCols; ++col) {
    Assert (isSubState(col) || isInsState(col) || isDelState(col), "Empty alignment column");
    if (inPath[col])
      ++inSeqIdx;
    if (outPath[col])
      ++outSeqIdx;
    col2InSeqIdx.push_back (inSeqIdx);
    col2OutSeqIdx.push_back (outSeqIdx);
  }
}

ForwardMatrix::ForwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : MutatorMatrix (mutatorParams, stock)
{
  sCell(0) = 0;
  for (AlignColIndex col = 1; col <= alignCols; ++col) {
    if (isSubState(col-1)) {
      sCell(col) = sCell(col-1) + mutatorScores.noGap + colSubScore(col-1);
    } else if (isInsState(col-1)) {
      for (Pos dupIdx = 0; dupIdx < maxDupLen - 1; ++dupIdx)
	tCell(col,dupIdx) = tCell(col-1,dupIdx+1) + colTanDupScore(col-1,dupIdx+1);
      sCell(col) = tCell(col-1,0) + colTanDupScore(col-1,0);
    } else if (isDelState(col-1)) {
      dCell(col) = log_sum_exp (sCell(col-1) + mutatorScores.delOpen,
				dCell(col-1) + mutatorScores.delExtend);
    }
    log_accum_exp (sCell(col), dCell(col) + mutatorScores.delEnd);
    for (Pos dupIdx = 0; dupIdx < maxDupLen - 1; ++dupIdx)
      log_accum_exp (tCell(col,dupIdx), sCell(col) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
  }
}

BackwardMatrix::BackwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : MutatorMatrix (mutatorParams, stock)
{
  sCell(alignCols) = 0;
  for (int col = alignCols; col >= 0; --col) {
    if (col < alignCols) {
      if (isSubState(col)) {
	sCell(col) = mutatorScores.noGap + colSubScore(col) + sCell(col+1);
      } else if (isInsState(col)) {
	for (Pos dupIdx = 1; dupIdx < maxDupLen; ++dupIdx)
	  tCell(col,dupIdx) = colTanDupScore(col,dupIdx) + tCell(col+1,dupIdx-1);
	tCell(col,0) = colTanDupScore(col,0) + sCell(col+1);
      } else if (isDelState(col)) {
	sCell(col) = mutatorScores.delOpen + dCell(col+1);
	dCell(col) = mutatorScores.delExtend + dCell(col+1);
      }
    }
    for (Pos dupIdx = 0; dupIdx < maxDupLen - 1; ++dupIdx)
      log_accum_exp (sCell(col), tCell(col,dupIdx) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
    log_accum_exp (dCell(col), sCell(col) + mutatorScores.delEnd);
  }
}

FwdBackMatrix::FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : fwd (mutatorParams, stock),
    back (mutatorParams, stock)
{ }

MutatorCounts FwdBackMatrix::counts() const {
  MutatorCounts counts (fwd.mutatorParams);
  const LogProb norm = fwd.loglike();
  for (AlignColIndex col = 1; col <= fwd.alignCols; ++col) {
    if (fwd.isSubState(col-1)) {
      const double c = exp (fwd.sCell(col-1) + fwd.mutatorScores.noGap + fwd.colSubScore(col-1) + back.sCell(col) - norm);
      counts.nNoGap += c;
      counts.nSub[fwd.colInBase(col)][fwd.colOutBase(col)] += c;
    } else if (fwd.isInsState(col-1)) {
      for (Pos dupIdx = 0; dupIdx < fwd.maxDupLen - 1; ++dupIdx) {
	const double ci = exp (fwd.tCell(col-1,dupIdx+1) + fwd.colTanDupScore(col-1,dupIdx+1) + back.tCell(col,dupIdx) - norm);
	counts.nSub[fwd.colTanDupBase(col-1,dupIdx+1)][fwd.colOutBase(col)] += ci;
      }
      const double c0 = exp (fwd.tCell(col-1,0) + fwd.colTanDupScore(col-1,0) + back.sCell(col) - norm);
      counts.nSub[fwd.colTanDupBase(col-1,0)][fwd.colOutBase(col)] += c0;
    } else if (fwd.isDelState(col-1)) {
      counts.nDelOpen += exp (fwd.sCell(col-1) + fwd.mutatorScores.delOpen + back.dCell(col) - norm);
      counts.nDelExtend += exp (fwd.dCell(col-1) + fwd.mutatorScores.delExtend + back.dCell(col) - norm);
    }
    counts.nDelEnd += exp (fwd.dCell(col) + fwd.mutatorScores.delEnd + back.sCell(col) - norm);
    for (Pos dupIdx = 0; dupIdx < fwd.maxDupLen - 1; ++dupIdx) {
      const double c = exp (fwd.sCell(col) + fwd.mutatorScores.tanDup + fwd.mutatorScores.len[dupIdx] + back.tCell(col,dupIdx));
      counts.nTanDup += c;
      counts.nLen[dupIdx] += c;
    }
  }
  return counts;
}

LogProb FwdBackMatrix::loglike() const {
  return fwd.loglike();
}

MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db) {
  MutatorParams current = init;
  LogProb best = -numeric_limits<double>::infinity();
  for (int iter = 0; iter < BaumWelchMaxIter; ++iter) {
    MutatorCounts counts (current);
    LogProb ll = 0;
    for (const auto& stock: db) {
      FwdBackMatrix fb (current, stock);
      counts += fb.counts();
      ll += fb.loglike();
    }
    ll += prior.logPrior (current);
    if ((ll - best) / abs(best) < BaumWelchMinFracInc)
      break;
    current = counts.mlParams (prior);
  }
  return current;
}
