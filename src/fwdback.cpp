#include <iomanip>
#include "fwdback.h"
#include "logsumexp.h"
#include "logger.h"

#define FwdBackTolerance 1e-5
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
  col2InSeqIdx.push_back (inSeqIdx);
  col2OutSeqIdx.push_back (outSeqIdx);
}

string MutatorMatrix::toString() const {
  ostringstream out;
  for (AlignColIndex col = 0; col <= alignCols; ++col) {
    out << setw(4) << col << ": "
	<< setw(10) << setprecision(6) << sCell(col) << "(S) "
	<< setw(10) << setprecision(6) << dCell(col) << "(D) ";
    for (Pos i = 0; i < maxDupLen; ++i)
      out << setw(10) << setprecision(6) << tCell(col,i) << "(T" << i+1 << ") ";
    out << "\n";
  }
  return out.str();
}

ForwardMatrix::ForwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : MutatorMatrix (mutatorParams, stock)
{
  sCell(0) = 0;
  for (AlignColIndex col = 1; col <= alignCols; ++col) {
    if (isSubState(col-1)) {
      sCell(col) = sCell(col-1) + mutatorScores.noGap + colSubScore(col-1);
    } else if (isInsState(col-1)) {
      for (Pos dupIdx = 0; dupIdx < maxDupLenAt(col-1) - 1; ++dupIdx)
	tCell(col,dupIdx) = tCell(col-1,dupIdx+1) + colTanDupScore(col-1,dupIdx+1);
      sCell(col) = tCell(col-1,0) + colTanDupScore(col-1,0);
    } else if (isDelState(col-1)) {
      dCell(col) = log_sum_exp (sCell(col-1) + mutatorScores.delOpen,
				dCell(col-1) + mutatorScores.delExtend);
    }
    log_accum_exp (sCell(col), dCell(col) + mutatorScores.delEnd);
    for (Pos dupIdx = 0; dupIdx < maxDupLenAt(col-1); ++dupIdx)
      log_accum_exp (tCell(col,dupIdx), sCell(col) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
  }
  LogThisAt(6,"Forward log-likelihood: " << loglike() << endl);
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
	for (Pos dupIdx = 1; dupIdx < maxDupLenAt(col); ++dupIdx)
	  tCell(col,dupIdx) = colTanDupScore(col,dupIdx) + tCell(col+1,dupIdx-1);
	tCell(col,0) = colTanDupScore(col,0) + sCell(col+1);
      } else if (isDelState(col)) {
	sCell(col) = mutatorScores.delOpen + dCell(col+1);
	dCell(col) = mutatorScores.delExtend + dCell(col+1);
      }
    }
    for (Pos dupIdx = 0; dupIdx < maxDupLenAt(col); ++dupIdx)
      log_accum_exp (sCell(col), tCell(col,dupIdx) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
    log_accum_exp (dCell(col), sCell(col) + mutatorScores.delEnd);
  }
  LogThisAt(6,"Backward log-likelihood: " << loglike() << endl);
}

FwdBackMatrix::FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock)
  : fwd (mutatorParams, stock),
    back (mutatorParams, stock)
{
  LogThisAt(9,"Forward matrix:\n" << fwd.toString() << "Backward matrix:\n" << back.toString());
  LogThisAt(8,"Forward-backward posterior probabilities:\n" << postProbsToString());
}

string FwdBackMatrix::postProbsToString() const {
  ostringstream out;
  for (AlignColIndex col = 1; col <= fwd.alignCols; ++col) {
    out << setw(4) << col << ": ";
    if (fwd.isSubState(col-1))
      out << setw(10) << setprecision(6) << pS2S(col) << "(S->S) ";
    else if (fwd.isInsState(col-1)) {
      for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(col-1) - 1; ++dupIdx)
	out << setw(10) << setprecision(6) << pT2T(col,dupIdx) << "(T" << dupIdx+2 << "->T" << dupIdx+1 << ") ";
      out << setw(10) << setprecision(6) << pT2S(col) << "(T1->S) ";
    } else if (fwd.isInsState(col-1))
      out << setw(10) << setprecision(6) << pS2D(col) << "(S->D) "
	  << setw(10) << setprecision(6) << pD2D(col) << "(D->D) ";
    for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(col-1); ++dupIdx)
      out << setw(10) << setprecision(6) << pS2T(col,dupIdx) << "(S->T" << dupIdx+1 << ") ";
    out << setw(10) << setprecision(6) << pD2S(col) << "(D->S)";
    out << "\n";
  }
  return out.str();
}

MutatorCounts FwdBackMatrix::counts() const {
  MutatorCounts counts (fwd.mutatorParams);
  for (AlignColIndex col = 1; col <= fwd.alignCols; ++col) {
    if (fwd.isSubState(col-1)) {
      const double c = pS2S(col);
      counts.nNoGap += c;
      counts.nSub[fwd.colInBase(col-1)][fwd.colOutBase(col-1)] += c;
    } else if (fwd.isInsState(col-1)) {
      for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(col-1) - 1; ++dupIdx) {
	const double ci = pT2T(col,dupIdx);
	counts.nSub[fwd.colTanDupBase(col-1,dupIdx+1)][fwd.colOutBase(col-1)] += ci;
      }
      const double c0 = pT2S(col);
      counts.nSub[fwd.colTanDupBase(col-1,0)][fwd.colOutBase(col-1)] += c0;
    } else if (fwd.isDelState(col-1)) {
      counts.nDelOpen += pS2D(col);
      counts.nDelExtend += pD2D(col);
    }
    counts.nDelEnd += pD2S(col);
    for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(col-1); ++dupIdx) {
      const double c = pS2T(col,dupIdx);
      counts.nTanDup += c;
      counts.nLen[dupIdx] += c;
    }
  }
  return counts;
}

LogProb FwdBackMatrix::loglike() const {
  if (abs ((fwd.loglike() - back.loglike()) / fwd.loglike()) > FwdBackTolerance)
    Warn ("Forward score (%g) does not match Backward score (%g)", fwd.loglike(), back.loglike());
  return fwd.loglike();
}

MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db) {
  MutatorParams current = init;
  LogProb best = -numeric_limits<double>::infinity();
  for (int iter = 0; iter < BaumWelchMaxIter; ++iter) {
    MutatorCounts counts (current);
    LogProb ll = 0;
    size_t nAlign = 0;
    for (const auto& stock: db) {
      FwdBackMatrix fb (current, stock);
      const auto stockCounts = fb.counts();
      const auto stockLoglike = fb.loglike();
      LogThisAt(5,"Counts for alignment #" << nAlign+1 << ":\n" << stockCounts.asJSON());
      LogThisAt(4,"Log-likelihood for alignment #" << nAlign+1 << ": " << stockLoglike << endl);
      counts += stockCounts;
      ll += stockLoglike;
      ++nAlign;
    }
    const LogProb lp = prior.logPrior (current);
    ll += lp;
    LogThisAt(6,"Log-prior: " << lp << endl);
    LogThisAt(2,"Iteration #" << iter+1 << ": log(likelihood*prior) = " << ll << endl);
    if ((ll - best) / abs(best) < BaumWelchMinFracInc)
      break;
    best = ll;
    LogThisAt(3,"Counts for iteration #" << iter+1 << ":\n" << counts.asJSON());
    current = counts.mlParams (prior);
    LogThisAt(5,"Parameters after iteration #" << iter+1 << ":\n" << current.asJSON());
  }
  return current;
}
