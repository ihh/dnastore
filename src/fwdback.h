#ifndef FWDBACK_INCLUDED
#define FWDBACK_INCLUDED

#include "mutator.h"
#include "stockholm.h"

class MutatorMatrix {
private:
  vguard<LogProb> cell;

  static inline size_t nCells (const MutatorParams& mutatorParams, const Stockholm& stock) {
    return (mutatorParams.maxDupLen() + 2) * (stock.columns() + 1);
  };

  inline size_t sCellIndex (AlignColIndex col) const {
    return (maxDupLen + 2) * col;
  };
  inline size_t dCellIndex (AlignColIndex col) const {
    return (maxDupLen + 2) * col + 1;
  };
  inline size_t tCellIndex (AlignColIndex col, Pos idx) const {
    return (maxDupLen + 2) * col + 2 + idx;
  };

protected:
  inline LogProb& sCell (AlignColIndex col) { return cell[sCellIndex(col)]; }
  inline LogProb& dCell (AlignColIndex col) { return cell[dCellIndex(col)]; }
  inline LogProb& tCell (AlignColIndex col, Pos idx) { return cell[tCellIndex(col,idx)]; }

public:
  const MutatorParams& mutatorParams;
  const MutatorScores mutatorScores;
  const size_t maxDupLen, alignCols;
  const Stockholm& stock;
  const Alignment align;
  const AlignRowPath& inPath;
  const AlignRowPath& outPath;
  TokSeq inSeq, outSeq;
  vguard<SeqIdx> col2InSeqIdx, col2OutSeqIdx;
  
  MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);

  inline LogProb sCell (AlignColIndex col) const { return cell[sCellIndex(col)]; }
  inline LogProb dCell (AlignColIndex col) const { return cell[dCellIndex(col)]; }
  inline LogProb tCell (AlignColIndex col, Pos idx) const { return cell[tCellIndex(col,idx)]; }

  inline bool isSubState (AlignColIndex col) const { return inPath[col] && outPath[col]; }
  inline bool isInsState (AlignColIndex col) const { return !inPath[col] && outPath[col]; }
  inline bool isDelState (AlignColIndex col) const { return inPath[col] && !outPath[col]; }

  inline Pos maxDupLenAt (AlignColIndex col) const { return min ((Pos) maxDupLen, (Pos) col2InSeqIdx[col]); }

  inline Base colInBase (AlignColIndex col) const { return inSeq[col2InSeqIdx[col]-1]; }
  inline Base colOutBase (AlignColIndex col) const { return outSeq[col2OutSeqIdx[col]-1]; }
  inline Base colTanDupBase (AlignColIndex col, Pos dupIdx) const { return inSeq[col2InSeqIdx[col]-1-dupIdx]; }

  inline LogProb colSubScore (AlignColIndex col) const {
    return mutatorScores.sub[colInBase(col)][colOutBase(col)];
  }

  inline LogProb colTanDupScore (AlignColIndex col, Pos dupIdx) const {
    return mutatorScores.sub[colTanDupBase(col,dupIdx)][colOutBase(col)];
  }

  string toString() const;
};

struct ForwardMatrix : MutatorMatrix {
  ForwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);
  inline LogProb loglike() const { return sCell(alignCols); }
};

struct BackwardMatrix : MutatorMatrix {
  BackwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);
  inline LogProb loglike() const { return sCell(0); }
};

struct FwdBackMatrix {
  ForwardMatrix fwd;
  BackwardMatrix back;
  FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);
  MutatorCounts counts() const;
  LogProb loglike() const;
  inline double pS2S (AlignColIndex destCol) const {
    return exp (fwd.sCell(destCol-1) + fwd.mutatorScores.noGap + fwd.colSubScore(destCol-1) + back.sCell(destCol) - loglike());
  }
  inline double pT2T (AlignColIndex destCol, Pos destDupIdx) const {
    return exp (fwd.tCell(destCol-1,destDupIdx+1) + fwd.colTanDupScore(destCol-1,destDupIdx+1) + back.tCell(destCol,destDupIdx) - loglike());
  }
  inline double pT2S (AlignColIndex destCol) const {
    return exp (fwd.tCell(destCol-1,0) + fwd.colTanDupScore(destCol-1,0) + back.sCell(destCol) - loglike());
  }
  inline double pS2D (AlignColIndex destCol) const {
    return exp (fwd.sCell(destCol-1) + fwd.mutatorScores.delOpen + back.dCell(destCol) - loglike());
  }
  inline double pD2D (AlignColIndex destCol) const {
    return exp (fwd.dCell(destCol-1) + fwd.mutatorScores.delExtend + back.dCell(destCol) - loglike());
  }
  inline double pD2S (AlignColIndex destCol) const {
    return exp (fwd.dCell(destCol) + fwd.mutatorScores.delEnd + back.sCell(destCol) - loglike());
  }
  inline double pS2T (AlignColIndex destCol, Pos destDupIdx) const {
    return exp (fwd.sCell(destCol) + fwd.mutatorScores.tanDup + fwd.mutatorScores.len[destDupIdx] + back.tCell(destCol,destDupIdx) - loglike());
  }
  string postProbsToString() const;
};

MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db);

#endif /* FWDBACK_INCLUDED */
