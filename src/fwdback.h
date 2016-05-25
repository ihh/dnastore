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
  inline size_t fCellIndex (AlignColIndex col, Pos idx) const {
    return (maxDupLen + 2) * col + 1 + idx;
  };

protected:
  const size_t maxDupLen, alignCols;

  inline LogProb& sCell (AlignColIndex col) { return cell[sCellIndex(col)]; }
  inline LogProb& dCell (AlignColIndex col) { return cell[dCellIndex(col)]; }
  inline LogProb& fCell (AlignColIndex col, Pos idx) { return cell[fCellIndex(col,idx)]; }

public:
  const MutatorParams& mutatorParams;
  const MutatorScores mutatorScores;

  const Stockholm& stock;
  const Alignment align;
  TokSeq inSeq, outSeq;
  vguard<SeqIdx> col2InSeqIdx, col2OutSeqIdx;
  
  MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);

  inline LogProb sCell (AlignColIndex col) const { return cell[sCellIndex(col)]; }
  inline LogProb dCell (AlignColIndex col) const { return cell[dCellIndex(col)]; }
  inline LogProb fCell (AlignColIndex col, Pos idx) const { return cell[fCellIndex(col,idx)]; }
};

struct FwdBackMatrix {
  MutatorMatrix fwd, back;
  FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock);
  MutatorCounts counts() const;
};

#endif /* FWDBACK_INCLUDED */
