#ifndef STOCKHOLM_INCLUDED
#define STOCKHOLM_INCLUDED

#include <map>
#include <list>
#include "fastseq.h"
#include "alignpath.h"

#define DefaultStockholmRowLength 80
#define MinStockholmCharsPerRow 10

struct Stockholm {
  vguard<FastSeq> gapped;
  map<string,string> gc;  // gc[tag][col]
  map<string,vguard<string> > gf;  // gf[tag][line]
  map<string,map<string,string> > gr;  // gr[tag][seqname][col]
  map<string,map<string,vguard<string> > > gs;  // gs[tag][seqname][line]

  Stockholm();
  Stockholm (istream& in);
  Stockholm (const vguard<FastSeq>& seq);
  
  void read (istream& in);
  void write (ostream& out, size_t charsPerRow = DefaultStockholmRowLength) const;

  size_t rows() const;
  size_t columns() const;
  AlignPath path() const;
};

list<Stockholm> readStockholmDatabase (const char* filename);

#endif /* STOCKHOLM_INCLUDED */
