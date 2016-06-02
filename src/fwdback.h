#ifndef FWDBACK_INCLUDED
#define FWDBACK_INCLUDED

#include "mutator.h"
#include "stockholm.h"

class MutatorMatrix {
private:
  vguard<LogProb> cell;

  static inline size_t nCells (const MutatorParams& mutatorParams, const Stockholm& stock) {
    const AlignPath path = stock.path();
    return (mutatorParams.maxDupLen() + 2) * (alignPathResiduesInRow (path.at(0)) + 1) * (alignPathResiduesInRow (path.at(1)) + 1);
  };

  inline size_t sCellIndex (SeqIdx inPos, SeqIdx outPos) const {
    return (maxDupLen + 2) * (inPos + (inLen + 1) * outPos);
  };
  inline size_t dCellIndex (SeqIdx inPos, SeqIdx outPos) const {
    return sCellIndex(inPos,outPos) + 1;
  };
  inline size_t tCellIndex (SeqIdx inPos, SeqIdx outPos, Pos dupIdx) const {
    return sCellIndex(inPos,outPos) + 2 + dupIdx;
  };

protected:
  inline LogProb& sCell (SeqIdx inPos, SeqIdx outPos) { return cell[sCellIndex(inPos,outPos)]; }
  inline LogProb& dCell (SeqIdx inPos, SeqIdx outPos) { return cell[dCellIndex(inPos,outPos)]; }
  inline LogProb& tCell (SeqIdx inPos, SeqIdx outPos, Pos idx) { return cell[tCellIndex(inPos,outPos,idx)]; }

public:
  const MutatorParams& mutatorParams;
  const MutatorScores mutatorScores;
  const size_t maxDupLen;
  const Stockholm& stock;
  const Alignment align;
  const GuideAlignmentEnvelope env;
  const TokSeq inSeq, outSeq;
  const size_t inLen, outLen;
  
  MutatorMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments);

  inline LogProb sCell (SeqIdx inPos, SeqIdx outPos) const { return cell[sCellIndex(inPos,outPos)]; }
  inline LogProb dCell (SeqIdx inPos, SeqIdx outPos) const { return cell[dCellIndex(inPos,outPos)]; }
  inline LogProb tCell (SeqIdx inPos, SeqIdx outPos, Pos dupIdx) const { return cell[tCellIndex(inPos,outPos,dupIdx)]; }

  inline Pos maxDupLenAt (SeqIdx inPos) const { return min ((Pos) maxDupLen, (Pos) inPos); }

  inline Base cellInBase (SeqIdx inPos) const { return inSeq[inPos-1]; }
  inline Base cellOutBase (SeqIdx outPos) const { return outSeq[outPos-1]; }
  inline Base cellTanDupBase (SeqIdx inPos, Pos dupIdx) const { return inSeq[inPos-1-dupIdx]; }

  inline LogProb cellSubScore (SeqIdx inPos, SeqIdx outPos) const {
    return mutatorScores.sub[cellInBase(inPos)][cellOutBase(outPos)];
  }

  inline LogProb cellTanDupScore (SeqIdx inPos, SeqIdx outPos, Pos dupIdx) const {
    return mutatorScores.sub[cellTanDupBase(inPos,dupIdx)][cellOutBase(outPos)];
  }

  string toString() const;
};

struct ForwardMatrix : MutatorMatrix {
  ForwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments);
  inline LogProb loglike() const { return sCell(inLen,outLen); }
};

struct BackwardMatrix : MutatorMatrix {
  BackwardMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments);
  inline LogProb loglike() const { return sCell(0,0); }
};

struct FwdBackMatrix {
  ForwardMatrix fwd;
  BackwardMatrix back;
  FwdBackMatrix (const MutatorParams& mutatorParams, const Stockholm& stock, bool strictAlignments);
  MutatorCounts counts() const;
  LogProb loglike() const;
  inline double pS2S (SeqIdx destInPos, SeqIdx destOutPos) const {
    return exp (fwd.sCell(destInPos-1,destOutPos-1) + fwd.mutatorScores.noGap + fwd.cellSubScore(destInPos,destOutPos) + back.sCell(destInPos,destOutPos) - loglike());
  }
  inline double pT2T (SeqIdx destInPos, SeqIdx destOutPos, Pos destDupIdx) const {
    return exp (fwd.tCell(destInPos,destOutPos-1,destDupIdx+1) + fwd.cellTanDupScore(destInPos,destOutPos,destDupIdx+1) + back.tCell(destInPos,destOutPos,destDupIdx) - loglike());
  }
  inline double pT2S (SeqIdx destInPos, SeqIdx destOutPos) const {
    return exp (fwd.tCell(destInPos,destOutPos-1,0) + fwd.cellTanDupScore(destInPos,destOutPos,0) + back.sCell(destInPos,destOutPos) - loglike());
  }
  inline double pS2D (SeqIdx destInPos, SeqIdx destOutPos) const {
    return exp (fwd.sCell(destInPos-1,destOutPos) + fwd.mutatorScores.delOpen + back.dCell(destInPos,destOutPos) - loglike());
  }
  inline double pD2D (SeqIdx destInPos, SeqIdx destOutPos) const {
    return exp (fwd.dCell(destInPos-1,destOutPos) + fwd.mutatorScores.delExtend + back.dCell(destInPos,destOutPos) - loglike());
  }
  inline double pD2S (SeqIdx destInPos, SeqIdx destOutPos) const {
    return exp (fwd.dCell(destInPos,destOutPos) + fwd.mutatorScores.delEnd + back.sCell(destInPos,destOutPos) - loglike());
  }
  inline double pS2T (SeqIdx destInPos, SeqIdx destOutPos, Pos destDupIdx) const {
    return exp (fwd.sCell(destInPos,destOutPos) + fwd.mutatorScores.tanDup + fwd.mutatorScores.len[destDupIdx] + back.tCell(destInPos,destOutPos,destDupIdx) - loglike());
  }
  string postProbsToString() const;
};

MutatorCounts expectedCounts (const MutatorParams& params, const list<Stockholm>& db, LogProb& ll, bool strictAlignments);
MutatorParams baumWelchParams (const MutatorParams& init, const MutatorCounts& prior, const list<Stockholm>& db, bool strictAlignments);

#endif /* FWDBACK_INCLUDED */
