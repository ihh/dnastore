#include <algorithm>
#include <iomanip>
#include "stockholm.h"
#include "regexmacros.h"
#include "util.h"

// POSIX basic regular expressions
const regex nonwhite_re (RE_DOT_STAR RE_NONWHITE_CHAR_CLASS RE_DOT_STAR, regex_constants::basic);
const regex seq_re (RE_WHITE_OR_EMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_OR_EMPTY, regex_constants::basic);
const regex gf_re (RE_WHITE_OR_EMPTY "#=GF" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_NONWHITE_CHAR_CLASS RE_DOT_STAR), regex_constants::basic);
const regex gc_re (RE_WHITE_OR_EMPTY "#=GC" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)), regex_constants::basic);
const regex gr_re (RE_WHITE_OR_EMPTY "#=GR" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)), regex_constants::basic);
const regex gs_re (RE_WHITE_OR_EMPTY "#=GS" RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_PLUS(RE_NONWHITE_CHAR_CLASS)) RE_WHITE_NONEMPTY RE_GROUP(RE_NONWHITE_CHAR_CLASS RE_DOT_STAR), regex_constants::basic);
const regex hash_re (RE_WHITE_OR_EMPTY "#" RE_DOT_STAR, regex_constants::basic);
const regex divider_re (RE_WHITE_OR_EMPTY "//" RE_WHITE_OR_EMPTY, regex_constants::basic);

Stockholm::Stockholm()
{ }

Stockholm::Stockholm (istream& in) {
  read (in);
}

Stockholm::Stockholm (const vguard<FastSeq>& seq)
  : gapped (seq)
{ }

void Stockholm::read (istream& in) {
  gf.clear();
  gc.clear();
  gs.clear();
  gr.clear();
  gapped.clear();

  smatch sm;
  map<string,string> seq;
  vguard<string> rowName;
  while (in && !in.eof()) {
    string line;
    getline(in,line);
    if (regex_match (line, sm, seq_re)) {
      if (!seq.count (sm.str(1)))
	rowName.push_back (sm.str(1));
      seq[sm.str(1)] += sm.str(2);
    } else if (regex_match (line, sm, gf_re))
      gf[sm.str(1)].push_back (sm.str(2));
    else if (regex_match (line, sm, gc_re))
      gc[sm.str(1)] += sm.str(2);
    else if (regex_match (line, sm, gr_re))
      gr[sm.str(2)][sm.str(1)] += sm.str(3);
    else if (regex_match (line, sm, gs_re))
      gs[sm.str(2)][sm.str(1)].push_back (sm.str(3));
    else if (regex_match (line, hash_re))
      continue;
    else if (regex_match (line, divider_re))
      break;
    else if (regex_match (line, nonwhite_re))
      Warn ("Unrecognized line in Stockholm file: %s", line.c_str());
  }
  for (const auto& name : rowName) {
    FastSeq fs;
    fs.name = name;
    fs.seq = seq[name];
    gapped.push_back (fs);
  }
}

void Stockholm::write (ostream& out, size_t charsPerRow) const {
  int nw = 0, tw = 0, w = 0, cols = columns();
  set<string> names;
  for (auto& fs : gapped) {
    w = max (w, (int) fs.name.size());
    names.insert (fs.name);
  }
  for (auto& tag_gf : gf)
    w = max (w, (int) tag_gf.first.size() + 5);
  for (auto& tag_gc : gc) {
    w = max (w, (int) tag_gc.first.size() + 5);
    cols = max (cols, (int) tag_gc.second.size());
  }
  for (auto& tag_gs : gs) {
    tw = max (tw, (int) tag_gs.first.size());
    for (auto& name_gs : tag_gs.second)
      nw = max (nw, (int) name_gs.first.size());
  }
  for (auto& tag_gr : gr) {
    tw = max (tw, (int) tag_gr.first.size());
    for (auto& name_gr : tag_gr.second) {
      nw = max (nw, (int) name_gr.first.size());
      cols = max (cols, (int) name_gr.second.size());
    }
  }
  if (tw > 0)
    w = max (w, nw + tw + 6);
  
  out << "# STOCKHOLM 1.0" << endl;
  for (auto& tag_gf : gf)
    for (auto& line : tag_gf.second)
      out << "#=GF " << left << setw(w-5) << tag_gf.first << " " << line << endl;

  for (auto& tag_gs : gs) {
    for (auto& fs : gapped)
      if (tag_gs.second.count (fs.name))
	for (auto& line : tag_gs.second.at(fs.name))
	  out << "#=GS " << left << setw(nw+1) << fs.name << left << setw(tw+1) << tag_gs.first << line << endl;
    for (auto& name_gs : tag_gs.second)
      if (!names.count (name_gs.first))
	for (auto& line : name_gs.second)
	  out << "#=GS " << left << setw(nw+1) << name_gs.first << left << setw(tw+1) << tag_gs.first << line << endl;
  }

  const int colStep = max (MinStockholmCharsPerRow, ((int) charsPerRow) - w - 1);
  for (int col = 0; col < cols; col += colStep) {
    for (auto& tag_gc : gc)
      if (col < tag_gc.second.size())
	out << "#=GC " << left << setw(w-5) << tag_gc.first << " " << tag_gc.second.substr(col,colStep) << endl;

    for (auto& fs : gapped) {
      if (col < fs.seq.size())
	out << left << setw(w+1) << fs.name << fs.seq.substr(col,colStep) << endl;
      for (auto& tag_gr : gr)
	if (tag_gr.second.count (fs.name))
	  if (col < tag_gr.second.at(fs.name).size())
	    out << "#=GR " << left << setw(nw+1) << fs.name << left << setw(tw+1) << tag_gr.first << tag_gr.second.at(fs.name).substr(col,colStep) << endl;
    }

    for (auto& tag_gr : gr)
      for (auto& name_gr : tag_gr.second)
	if (!names.count (name_gr.first))
	  if (col < name_gr.second.size())
	    out << "#=GR " << left << setw(nw+1) << name_gr.first << left << setw(tw+1) << tag_gr.first << name_gr.second.substr(col,colStep) << endl;

    if (col + colStep < cols)
      out << endl;
  }
  out << "//" << endl;
}

size_t Stockholm::rows() const {
  return gapped.size();
}

size_t Stockholm::columns() const {
  return gappedSeqColumns (gapped);
}

AlignPath Stockholm::path() const {
  Alignment a (gapped);
  return a.path;
}
