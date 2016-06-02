#include <iomanip>
#include "fwdback.h"
#include "logsumexp.h"
#include "logger.h"

#define FwdBackTolerance 1e-5
#define BaumWelchMinFracInc .001
#define BaumWelchMaxIter 100

MutatorMatrix::MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments)
  : cell (nCells (mutatorParams, stock), -numeric_limits<double>::infinity()),
    mutatorParams (mutatorParams),
    mutatorScores (mutatorParams),
    maxDupLen (mutatorParams.maxDupLen()),
    stock (stock),
    align (stock.gapped),
    env (align.path, 0, 1, strictAlignments ? 0 : 2*maxDupLen),
    inSeq (align.ungapped.at(0).tokens(dnaAlphabetString)),
    outSeq (align.ungapped.at(1).tokens(dnaAlphabetString)),
    inLen (inSeq.size()),
    outLen (outSeq.size())
{
  Assert (stock.rows() == 2, "Training mutator model requires a 2-row alignment; this alignment has %d rows", stock.rows());
}

string MutatorMatrix::toString() const {
  ostringstream out;
  for (SeqIdx ip = 0; ip <= inLen; ++ip)
    for (SeqIdx op = 0; op <= outLen; ++op)
      if (env.inRange(ip,op)) {
	out << setw(4) << ip << setw(4) << op << ": "
	    << setw(10) << setprecision(5) << sCell(ip,op) << "(S) "
	    << setw(10) << setprecision(5) << dCell(ip,op) << "(D) ";
	for (Pos i = 0; i < maxDupLen; ++i)
	  out << setw(10) << setprecision(5) << tCell(ip,op,i) << "(T" << i+1 << ") ";
	out << "\n";
      }
  return out.str();
}

ForwardMatrix::ForwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments)
  : MutatorMatrix (mutatorParams, stock, strictAlignments)
{
  sCell(0,0) = 0;
  for (SeqIdx ip = 0; ip <= inLen; ++ip)
    for (SeqIdx op = 0; op <= outLen; ++op)
      if (env.inRange(ip,op)) {
	if (ip > 0 && op > 0 && env.inRange(ip-1,op-1))
	  sCell(ip,op) = sCell(ip-1,op-1) + mutatorScores.noGap + cellSubScore(ip,op);
	if (ip > 0 && op > 0 && env.inRange(ip,op-1)) {
	  for (Pos dupIdx = 0; dupIdx < maxDupLenAt(ip) - 1; ++dupIdx)
	    tCell(ip,op,dupIdx) = tCell(ip,op-1,dupIdx+1) + cellTanDupScore(ip,op,dupIdx+1);
	  log_accum_exp (sCell(ip,op), tCell(ip,op-1,0) + cellTanDupScore(ip,op,0));
	}
	if (ip > 0 && env.inRange(ip-1,op))
	  dCell(ip,op) = log_sum_exp (sCell(ip-1,op) + mutatorScores.delOpen,
				      dCell(ip-1,op) + mutatorScores.delExtend);
	log_accum_exp (sCell(ip,op), dCell(ip,op) + mutatorScores.delEnd);
	for (Pos dupIdx = 0; dupIdx < maxDupLenAt(ip); ++dupIdx)
	  log_accum_exp (tCell(ip,op,dupIdx), sCell(ip,op) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
      }
  LogThisAt(6,"Forward log-odds ratio: " << loglike() << endl);
}

BackwardMatrix::BackwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments)
  : MutatorMatrix (mutatorParams, stock, strictAlignments)
{
  sCell(inLen,outLen) = 0;
  for (int ip = inLen; ip >= 0; --ip)
    for (int op = outLen; op >= 0; --op)
      if (env.inRange(ip,op)) {
	if (ip < inLen && op < outLen && env.inRange(ip+1,op+1))
	  sCell(ip,op) = mutatorScores.noGap + cellSubScore(ip+1,op+1) + sCell(ip+1,op+1);
	if (ip > 0 && op < outLen && env.inRange(ip,op+1)) {
	  for (Pos dupIdx = 1; dupIdx < maxDupLenAt(ip); ++dupIdx)
	    tCell(ip,op,dupIdx) = cellTanDupScore(ip,op+1,dupIdx) + tCell(ip,op+1,dupIdx-1);
	  tCell(ip,op,0) = cellTanDupScore(ip,op+1,0) + sCell(ip,op+1);
	}
	if (ip < inLen && env.inRange(ip+1,op)) {
	  log_accum_exp (sCell(ip,op), mutatorScores.delOpen + dCell(ip+1,op));
	  dCell(ip,op) = mutatorScores.delExtend + dCell(ip+1,op);
	}
	for (Pos dupIdx = 0; dupIdx < maxDupLenAt(ip); ++dupIdx)
	  log_accum_exp (sCell(ip,op), tCell(ip,op,dupIdx) + mutatorScores.tanDup + mutatorScores.len[dupIdx]);
	log_accum_exp (dCell(ip,op), sCell(ip,op) + mutatorScores.delEnd);
      }
  LogThisAt(6,"Backward log-odds ratio: " << loglike() << endl);
}

FwdBackMatrix::FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments)
  : fwd (mutatorParams, stock, strictAlignments),
    back (mutatorParams, stock, strictAlignments)
{
  LogThisAt(7,"Scores:\n" << fwd.mutatorScores.toJSON());
  LogThisAt(9,"Forward matrix:\n" << fwd.toString() << "Backward matrix:\n" << back.toString());
  LogThisAt(8,"Forward-backward posterior probabilities:\n" << postProbsToString());

  if (abs ((fwd.loglike() - back.loglike()) / fwd.loglike()) > FwdBackTolerance)
    Warn ("Forward score (%g) does not match Backward score (%g)", fwd.loglike(), back.loglike());
}

string FwdBackMatrix::postProbsToString() const {
  ostringstream out;
  for (SeqIdx ip = 0; ip <= fwd.inLen; ++ip)
    for (SeqIdx op = 0; op <= fwd.outLen; ++op)
      if (fwd.env.inRange(ip,op)) {
	out << setw(4) << ip << setw(4) << op << ": ";
	if (ip > 0 && op > 0)
	  out << setw(10) << setprecision(5) << pS2S(ip,op) << "(S->S) ";
	if (ip > 0)
	  out << setw(10) << setprecision(5) << pS2D(ip,op) << "(S->D) "
	      << setw(10) << setprecision(5) << pD2D(ip,op) << "(D->D) ";
	if (ip > 0 && op > 0) {
	  out << setw(10) << setprecision(5) << pT2S(ip,op) << "(T1->S) ";
	  for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(ip) - 1; ++dupIdx)
	    out << setw(10) << setprecision(5) << pT2T(ip,op,dupIdx) << "(T" << dupIdx+2 << "->T" << dupIdx+1 << ") ";
	}
	for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(ip); ++dupIdx)
	  out << setw(10) << setprecision(5) << pS2T(ip,op,dupIdx) << "(S->T" << dupIdx+1 << ") ";
	out << setw(10) << setprecision(5) << pD2S(ip,op) << "(D->S)";
	out << "\n";
    }
  return out.str();
}

MutatorCounts FwdBackMatrix::counts() const {
  MutatorCounts counts (fwd.mutatorParams);
  for (SeqIdx ip = 0; ip <= fwd.inLen; ++ip)
    for (SeqIdx op = 0; op <= fwd.outLen; ++op)
      if (fwd.env.inRange(ip,op)) {
	if (ip > 0 && op > 0) {
	  const double c = pS2S(ip,op);
	  counts.nNoGap += c;
	  counts.nSub[fwd.cellInBase(ip)][fwd.cellOutBase(op)] += c;
	}
	if (ip > 0 && op > 0) {
	  for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(ip) - 1; ++dupIdx) {
	    const double ci = pT2T(ip,op,dupIdx);
	    counts.nSub[fwd.cellTanDupBase(ip,dupIdx+1)][fwd.cellOutBase(op)] += ci;
	  }
	  const double c0 = pT2S(ip,op);
	  counts.nSub[fwd.cellTanDupBase(ip,0)][fwd.cellOutBase(op)] += c0;
	}
	if (ip > 0) {
	  counts.nDelOpen += pS2D(ip,op);
	  counts.nDelExtend += pD2D(ip,op);
	}
	counts.nDelEnd += pD2S(ip,op);
	for (Pos dupIdx = 0; dupIdx < fwd.maxDupLenAt(ip); ++dupIdx) {
	  const double c = pS2T(ip,op,dupIdx);
	  counts.nTanDup += c;
	  counts.nLen[dupIdx] += c;
	}
      }
  return counts;
}

LogProb FwdBackMatrix::loglike() const {
  return fwd.loglike();
}

MutatorCounts expectedCounts (const MutatorParams& params, const list<Stockholm>& db, LogProb& ll, bool strictAlignments) {
  MutatorCounts counts (params);
  ll = 0;
  size_t nAlign = 0;
  const size_t nTotal = db.size();
  ProgressLog (plog, 2);
  plog.initProgress ("Getting Baum-Welch counts (%u alignments)", nTotal);
  for (const auto& stock: db) {
    plog.logProgress (nAlign / (double) nTotal, "sequence %u/%u", nAlign+1, nTotal);
    FwdBackMatrix fb (params, stock, strictAlignments);
    const auto stockCounts = fb.counts();
    const auto stockLoglike = fb.loglike();
    LogThisAt(5,"Counts for alignment #" << nAlign+1 << ":\n" << stockCounts.asJSON());
    LogThisAt(4,"Log-odds ratio for alignment #" << nAlign+1 << ": " << stockLoglike << endl);
    counts += stockCounts;
    ll += stockLoglike;
    ++nAlign;
  }
  return counts;
}

MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db, bool strictAlignments) {
  MutatorParams current = init;
  LogProb best = -numeric_limits<double>::infinity();
  for (int iter = 0; iter < BaumWelchMaxIter; ++iter) {
    LogProb ll;
    const MutatorCounts counts = expectedCounts (current, db, ll, strictAlignments);
    const LogProb lp = prior.logPrior (current);
    ll += lp;
    LogThisAt(6,"Log-prior: " << lp << endl);
    LogThisAt(2,"Iteration #" << iter+1 << ": log(oddsRatio*prior) = " << ll << endl);
    if ((ll - best) / abs(best) < BaumWelchMinFracInc)
      break;
    best = ll;
    LogThisAt(3,"Counts for iteration #" << iter+1 << ":\n" << counts.asJSON());
    current = counts.mlParams (prior);
    current.local = init.local;
    LogThisAt(5,"Parameters after iteration #" << iter+1 << ":\n" << current.asJSON());
  }
  return current;
}
