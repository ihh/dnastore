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
    return (maxDupLen + 2) * col + 1 + idx;
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

  inline Base colInBase (AlignColIndex col) const { return inSeq[col2InSeqIdx[col]-1]; }
  inline Base colOutBase (AlignColIndex col) const { return outSeq[col2OutSeqIdx[col]-1]; }
  inline Base colTanDupBase (AlignColIndex col, Pos dupIdx) const { return inSeq[col2InSeqIdx[col]-1-dupIdx]; }

  inline LogProb colSubScore (AlignColIndex col) const {
    return mutatorScores.sub[colInBase(col)][colOutBase(col)];
  }

  inline LogProb colTanDupScore (AlignColIndex col, Pos dupIdx) const {
    return mutatorScores.sub[colTanDupBase(col,dupIdx)][colOutBase(col)];
  }
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
};

MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db);

#endif /* FWDBACK_INCLUDED */
