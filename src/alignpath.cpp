#include "alignpath.h"
#include "util.h"
#include "logger.h"

const char Alignment::gapChar = '-';
const char Alignment::wildcardChar = '*';

// map used by alignPathMerge
struct AlignSeqMap {
  typedef size_t AlignNum;
  const vguard<AlignPath>& alignments;
  map<AlignRowIndex,SeqIdx> seqLen;
  vguard<AlignColIndex> alignCols;
  map<AlignNum,map<AlignColIndex,map<AlignRowIndex,SeqIdx> > > alignColRowToPos;
  map<AlignRowIndex,map<SeqIdx,map<AlignNum,AlignColIndex> > > rowPosAlignToCol;
  AlignSeqMap (const vguard<AlignPath>& alignments);
  map<AlignNum,AlignColIndex> linkedColumns (AlignNum nAlign, AlignColIndex col) const;
};

AlignColIndex gappedSeqColumns (const vguard<FastSeq>& gapped) {
  AlignColIndex cols = 0;
  for (size_t row = 0; row < gapped.size(); ++row)
    if (row == 0)
      cols = gapped[row].length();
    else
      Assert (cols == gapped[row].length(), "Alignment is not flush: sequence %s has %u chars, but sequence %s has %u chars", gapped[0].name.c_str(), cols, gapped[row].name.c_str(), gapped[row].length());

  return cols;
}

AlignColIndex alignPathColumns (const AlignPath& a) {
  AlignColIndex cols = 0;
  bool first = true;
  AlignRowIndex firstRow = 0;
  for (auto& row_path : a) {
    if (first) {
      firstRow = row_path.first;
      cols = row_path.second.size();
      first = false;
    } else
      Assert (cols == row_path.second.size(), "Alignment path is not flush: row %u has %u columns, but row %u has %u columns", firstRow, cols, row_path.first, row_path.second.size());
  }

  return cols;
}

SeqIdx alignPathResiduesInRow (const AlignRowPath& r) {
  SeqIdx l = 0;
  for (bool b : r)
    if (b)
      ++l;
  return l;
}

AlignPath alignPathUnion (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  a.insert (a2.begin(), a2.end());
  return a;
}

AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2) {
  AlignPath a = a1;
  const AlignColIndex c1 = alignPathColumns(a1), c2 = alignPathColumns(a2);
  for (auto& iter : a)
    if (a2.find(iter.first) == a2.end())
      iter.second.insert (iter.second.end(), c2, false);
  for (auto& iter2 : a2) {
    const AlignRowIndex row = iter2.first;
    const AlignRowPath& rPath = iter2.second;
    AlignRowPath& lPath = a[row];
    if (lPath.empty())
      lPath.insert (lPath.end(), c1, false);
    lPath.insert (lPath.end(), rPath.begin(), rPath.end());
  }
  return a;
}

AlignPath alignPathConcat (const AlignPath& a1, const AlignPath& a2, const AlignPath& a3) {
  return alignPathConcat (alignPathConcat (a1, a2), a3);
}

AlignSeqMap::AlignSeqMap (const vguard<AlignPath>& alignments)
  : alignments (alignments)
{
  // get row indices and sequence lengths; confirm row & sequence lengths match
  for (auto& align : alignments) {
    if (align.size() == 0)
      alignCols.push_back (0);
    else {
      alignCols.push_back (alignPathColumns (align));
      for (auto& row_path : align) {
	const AlignRowIndex row = row_path.first;
	const AlignRowPath& path = row_path.second;
	SeqIdx len = alignPathResiduesInRow (path);
	if (seqLen.find(row) == seqLen.end())
	  seqLen[row] = len;
	else
	  Assert (seqLen[row] == len, "Incompatible number of residues for row #%d of alignment (%d != %d)", row, seqLen[row], len);
      }
    }
  }

  // build bidirectional map from (align#,column#) <==> (row#,residue#)
  for (size_t nAlign = 0; nAlign < alignments.size(); ++nAlign) {
    auto& align = alignments[nAlign];
    map<AlignRowIndex,SeqIdx> rowPos;
    for (auto& row_path : align)
      rowPos[row_path.first] = 0;
    for (AlignColIndex col = 0; col < alignCols[nAlign]; ++col) {
      for (auto& row_path : align)
	if (row_path.second[col]) {
	  const SeqIdx pos = rowPos[row_path.first]++;
	  alignColRowToPos[nAlign][col][row_path.first] = pos;
	  rowPosAlignToCol[row_path.first][pos][nAlign] = col;
	}
    }
  }
}

map<AlignSeqMap::AlignNum,AlignColIndex> AlignSeqMap::linkedColumns (AlignNum nAlign, AlignColIndex col) const {
  map<AlignNum,AlignColIndex> ac, acQueue;
  acQueue[nAlign] = col;
  while (acQueue.size() > ac.size()) {
    for (auto& nAlign_col : acQueue)
      if (ac.find(nAlign_col.first) == ac.end()) {
	ac.insert (nAlign_col);
	for (auto& row_pos : alignColRowToPos.at(nAlign_col.first).at(nAlign_col.second))
	  for (auto& linked_nAlign_col : rowPosAlignToCol.at(row_pos.first).at(row_pos.second)) {
	    if (ac.find (linked_nAlign_col.first) != ac.end())
	      Assert (ac[linked_nAlign_col.first] == linked_nAlign_col.second, "Inconsistent alignments\nColumn %u of alignment %u points to position %u of sequence %u, which points back to column %u of alignment %u", col, nAlign, row_pos.second, row_pos.first, linked_nAlign_col.second, linked_nAlign_col.first);
	    acQueue.insert (linked_nAlign_col);
	  }
      }
  }
  return ac;
}

AlignPath alignPathMerge (const vguard<AlignPath>& alignments) {
  const AlignSeqMap alignSeqMap (alignments);
  AlignPath a;
  for (auto& row_seqlen : alignSeqMap.seqLen)
    a[row_seqlen.first].clear();
  vguard<AlignColIndex> nextCol (alignments.size(), 0);
  bool allDone, noneReady;
  do {
    allDone = noneReady = true;
    map<AlignSeqMap::AlignNum,AlignColIndex> linkedCols;
    for (AlignSeqMap::AlignNum n = 0; n < alignments.size(); ++n)
      if (nextCol[n] < alignSeqMap.alignCols[n]) {
	allDone = false;
	bool ready = true;
	linkedCols = alignSeqMap.linkedColumns (n, nextCol[n]);
	for (const auto& nAlign_col : linkedCols) {
	  if (nextCol[nAlign_col.first] != nAlign_col.second) {
	    ready = false;
	    break;
	  }
	}
	if (ready) {
	  noneReady = false;
	  if (linkedCols.size()) {
	    for (auto& idx_path : a)
	      idx_path.second.push_back (false);
	    for (const auto& nAlign_col : linkedCols) {
	      for (const auto& row_path : alignments.at(nAlign_col.first))
		if (alignments.at(nAlign_col.first).at(row_path.first).at(nAlign_col.second))
		  a[row_path.first].back() = true;
	      ++nextCol[nAlign_col.first];
	    }
	  } else
	    ++nextCol[n];  // empty column
	  break;
	}
      }
    if (noneReady && !allDone) {
      for (AlignSeqMap::AlignNum n = 0; n < alignments.size(); ++n)
	cerr << "Alignment #" << n << ": next column " << nextCol[n] << endl;
      Abort ("%s fail, no alignments ready", __func__);
    }
  } while (!allDone);

  const AlignRowIndex rows = a.size();
  const AlignColIndex cols = alignPathColumns (a);  // this will also test if alignment is flush
  LogThisAt(2,"Merged " << alignments.size() << " alignments into a single alignment with " << rows << " rows and " << cols << " columns" << endl);

  return a;
}

Alignment::Alignment (const vguard<FastSeq>& gapped)
  : ungapped (gapped.size())
 {
  for (AlignRowIndex row = 0; row < gapped.size(); ++row) {
    ungapped[row].name = gapped[row].name;
    ungapped[row].comment = gapped[row].comment;
    AlignRowPath rowPath (gapped[row].length(), false);
    for (AlignColIndex col = 0; col < rowPath.size(); ++col)
      if (!isGap (gapped[row].seq[col])) {
	rowPath[col] = true;
	ungapped[row].seq.push_back (gapped[row].seq[col]);
	ungapped[row].qual.push_back (gapped[row].qual[col]);
      }
    path[row] = rowPath;
  }
}

Alignment::Alignment (const vguard<FastSeq>& ungapped, const AlignPath& path)
  : ungapped(ungapped), path(path)
{ }

vguard<FastSeq> Alignment::gapped() const {
  vguard<FastSeq> gs (ungapped.size());
  for (auto& row_path : path) {
    FastSeq& g = gs[row_path.first];
    const FastSeq& ug = ungapped[row_path.first];
    const AlignColIndex cols = row_path.second.size();
    g.name = ug.name;
    g.comment = ug.comment;
    g.seq.reserve (cols);
    g.qual.reserve (cols);
    SeqIdx pos = 0;
    for (AlignColIndex col = 0; col < cols; ++col)
      if (row_path.second[col]) {
	Assert (ug.seq.size() > pos, "Sequence position %u out of bounds for sequence %s", col, ug.name.c_str());
	g.seq.push_back (ug.seq[pos]);
	if (ug.hasQual()) {
	  Assert (ug.qual.size() > pos, "Quality score at position %u out of bounds for sequence %s", col, ug.name.c_str());
	  g.qual.push_back (ug.qual[pos]);
	}
	++pos;
      } else {
	g.seq.push_back ('-');
	g.qual.push_back ('!');
      }
  }
  return gs;
}
GuideAlignmentEnvelope::GuideAlignmentEnvelope (const AlignPath& guide, AlignRowIndex row1, AlignRowIndex row2, int maxDistance)
  : maxDistance (maxDistance),
    row1 (row1),
    row2 (row2)
{
  Assert (guide.find(row1) != guide.end(), "Guide alignment is missing row #%u", row1);
  Assert (guide.find(row2) != guide.end(), "Guide alignment is missing row #%u", row2);

  const AlignColIndex cols = alignPathColumns (guide);
  cumulativeMatches.reserve (cols + 1);
  int matches = 0;

  row1PosToCol.push_back (0);
  row2PosToCol.push_back (0);
  cumulativeMatches.push_back (0);
  
  for (AlignColIndex col = 0; col < cols; ++col) {
    if (guide.at(row1)[col])
      row1PosToCol.push_back (col + 1);

    if (guide.at(row2)[col])
      row2PosToCol.push_back (col + 1);

    if (guide.at(row1)[col] && guide.at(row2)[col])
      ++matches;

    cumulativeMatches.push_back (matches);
  }
}
