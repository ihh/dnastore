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

#endif /* ALIGNPATH_INCLUDED */
