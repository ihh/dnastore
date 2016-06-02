#ifndef ALIGNPATH_INCLUDED
#define ALIGNPATH_INCLUDED

#include <map>
#include <set>
#include "vguard.h"
#include "fastseq.h"

typedef size_t AlignRowIndex;
typedef size_t AlignColIndex;
typedef vguard<bool> AlignRowPath;
typedef map<AlignRowIndex,AlignRowPath> AlignPath;

AlignColIndex gappedSeqColumns (const vguard<FastSeq>& gapped);
AlignColIndex alignPathColumns (const AlignPath& a);
SeqIdx alignPathResiduesInRow (const AlignRowPath& r);

AlignPath alignPathUnion (const AlignPath& a1, const AlignPath& a2);  // simple union (no AlignRowIndex shared between a1 & a2)
AlignPath alignPathMerge (const vguard<AlignPath>& alignments);  // synchronized merge
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2);  // lengthwise concatenation
AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2, const AlignPath& a3);

struct Alignment {
  static const char gapChar, wildcardChar;
  vguard<FastSeq> ungapped;
  AlignPath path;
  Alignment() { }
  Alignment (const vguard<FastSeq>& gapped);
  Alignment (const vguard<FastSeq>& ungapped, const AlignPath& path);
  vguard<FastSeq> gapped() const;
  static inline bool isGap (char c) { return c == '-' || c == '.'; }
  static inline bool isWildcard (char c) { return c == wildcardChar; }
};

struct GuideAlignmentEnvelope {
  // cumulativeMatches[col] = number of matches before column #col of pairwise alignment of (row1,row2)
  vguard<int> cumulativeMatches;
  // rowPosToCol[0 or 1][seqpos] = alignment column number of position #seqpos of (row1 or row2)
  vguard<AlignColIndex> row1PosToCol, row2PosToCol;
  AlignRowIndex row1, row2;
  int maxDistance;

  GuideAlignmentEnvelope() : maxDistance(-1) { }
  GuideAlignmentEnvelope (const AlignPath& guide, AlignRowIndex row1, AlignRowIndex row2, int maxDistance);

  inline bool initialized() const { return maxDistance >= 0; }

  inline bool inRange (SeqIdx pos1, SeqIdx pos2) const {
    if (!initialized())
      return true;
    const int d = cumulativeMatches[row1PosToCol[pos1]] - cumulativeMatches[row2PosToCol[pos2]];
    return abs(d) <= maxDistance;
  }
};

#endif /* ALIGNPATH_INCLUDED */
