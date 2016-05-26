#include <fstream>
#include "mutator.h"
#include "jsonutil.h"
#include "logger.h"

void MutatorParams::writeJSON (ostream& out) const {
  out << "{\n";
  out << " \"pDelOpen\": " << pDelOpen << ",\n";
  out << " \"pDelExtend\": " << pDelExtend << ",\n";
  out << " \"pTanDup\": " << pTanDup << ",\n";
  out << " \"pTransition\": " << pTransition << ",\n";
  out << " \"pTransversion\": " << pTransversion << ",\n";
  out << " \"pLen\": [ " << to_string_join(pLen,", ") << " ],\n";
  out << " \"local\": " << (local ? "true" : "false") << "\n";
  out << "}\n";
}

void MutatorParams::readJSON (istream& in) {
  pLen.clear();
  ParsedJson pj (in);
  pDelOpen = pj.getNumber ("pDelOpen");
  pDelExtend = pj.getNumber ("pDelExtend");
  pTanDup = pj.getNumber ("pTanDup");
  pTransition = pj.getNumber ("pTransition");
  pTransversion = pj.getNumber ("pTransversion");
  local = pj.getBool("local");
  JsonValue pLenArray = pj.getType ("pLen", JSON_ARRAY);
  for (JsonIterator iter = begin(pLenArray); iter != end(pLenArray); ++iter)
    pLen.push_back (iter->value.toNumber());
}

string MutatorParams::asJSON() const {
  ostringstream out;
  writeJSON (out);
  return out.str();
}

MutatorParams MutatorParams::fromJSON (istream& in) {
  MutatorParams mp;
  mp.readJSON (in);
  return mp;
}

MutatorParams MutatorParams::fromFile (const char* filename) {
  ifstream infile (filename);
  if (!infile)
    Fail ("File not found: %s", filename);
  return fromJSON (infile);
}

MutatorParams& MutatorParams::initMaxDupLen (size_t maxDupLen) {
  pLen = vguard<double> (maxDupLen, 1. / (double) maxDupLen);
  return *this;
}

MutatorScores::MutatorScores (const MutatorParams& params)
  : delOpen (log (params.pDelOpen)),
    tanDup (log (params.pTanDup)),
    noGap (log (params.pNoGap())),
    delExtend (log (params.pDelExtend)),
    delEnd (log (params.pDelEnd())),
    sub (4, vguard<LogProb> (4)),
    len (params.maxDupLen())
{
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      sub[i][j] = i==j
	? log(params.pMatch())
	: (isTransition(i,j)
	   ? log(params.pTransition)
	   : log(params.pTransversion/2));
  for (Pos l = 0; l < params.maxDupLen(); ++l)
    len[l] = log(params.pLen[l]);
}

MutatorCounts::MutatorCounts (const MutatorParams& params)
  : nDelOpen(0),
    nTanDup(0),
    nNoGap(0),
    nDelExtend(0),
    nDelEnd(0),
    nSub(4,vguard<double>(4,0)),
    nLen(params.maxDupLen(),0)
{ }

void MutatorCounts::writeJSON (ostream& out) const {
  out << "{\n";
  out << " \"nDelOpen\": " << nDelOpen << ",\n";
  out << " \"nTanDup\": " << nTanDup << ",\n";
  out << " \"nNoGap\": " << nNoGap << ",\n";
  out << " \"nDelExtend\": " << nDelExtend << ",\n";
  out << " \"nDelEnd\": " << nDelEnd << ",\n";
  out << " \"nLen\": [ " << to_string_join(nLen,", ") << " ],\n";
  out << " \"nSub\": [ ";
  for (Base i = 0; i < 4; ++i)
    out << (i > 0 ? ", " : "") << "[" << to_string_join(nSub[i],",") << "]";
  out << " ],\n";
  out << " \"nMatch\": " << nMatch() << ",\n";
  out << " \"nTransition\": " << nTransition() << ",\n";
  out << " \"nTransversion\": " << nTransversion() << "\n";
  out << "}\n";
}

string MutatorCounts::asJSON() const {
  ostringstream out;
  writeJSON (out);
  return out.str();
}

MutatorCounts& MutatorCounts::initLaplace (double n) {
  nDelOpen = n;
  nTanDup = n;
  nNoGap = n;
  nDelExtend = n;
  nDelEnd = n;
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      nSub[i][j] = n;
  for (size_t l = 0; l < nLen.size(); ++l)
    nLen[l] = n;
  return *this;
}

MutatorCounts& MutatorCounts::operator+= (const MutatorCounts& c) {
  Assert (nLen.size() == c.nLen.size(), "Length mismatch");
  nDelOpen += c.nDelOpen;
  nTanDup += c.nTanDup;
  nNoGap += c.nNoGap;
  nDelExtend += c.nDelExtend;
  nDelEnd += c.nDelEnd;
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      nSub[i][j] += c.nSub[i][j];
  for (size_t l = 0; l < nLen.size(); ++l)
    nLen[l] += c.nLen[l];
  return *this;
}

MutatorCounts MutatorCounts::operator+ (const MutatorCounts& c) const {
  MutatorCounts r (*this);
  r += c;
  return r;
}

MutatorParams MutatorCounts::mlParams() const {
  MutatorParams p;
  p.initMaxDupLen (nLen.size());
  p.pDelOpen = nDelOpen / (nDelOpen + nTanDup + nNoGap);
  p.pTanDup = nTanDup / (nDelOpen + nTanDup + nNoGap);
  p.pDelExtend = nDelExtend / (nDelExtend + nDelEnd);
  const double ni = nTransition(), nv = nTransversion(), nm = nMatch();
  p.pTransition = ni / (ni + nv + nm);
  p.pTransversion = nv / (ni + nv + nm);
  return p;
}

double MutatorCounts::nMatch() const {
  double n = 0;
  for (Base i = 0; i < 4; ++i)
    n += nSub[i][i];
  return n;
}

double MutatorCounts::nTransition() const {
  double n = 0;
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      if (isTransition(i,j))
	n += nSub[i][j];
  return n;
}

double MutatorCounts::nTransversion() const {
  double n = 0;
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      if (isTransversion(i,j))
	n += nSub[i][j];
  return n;
}

MutatorParams MutatorCounts::mlParams (const MutatorCounts& prior) const {
  const MutatorCounts c = *this + prior;
  LogThisAt(10,"Counts + prior:\n" << c.asJSON());
  return c.mlParams();
}

LogProb MutatorCounts::logPrior (const MutatorParams& params) const {
  const vguard<double> pGap = { params.pDelOpen, params.pTanDup, params.pNoGap() };
  const vguard<double> nGap = { nDelOpen, nTanDup, nNoGap };

  const vguard<double> pSub = { params.pTransition, params.pTransversion, params.pMatch() };
  const vguard<double> nSub = { nTransition(), nTransversion(), nMatch() };

  return logBetaPdfCounts(params.pDelExtend,nDelExtend,nDelEnd)
    + logDirichletPdfCounts(pGap,nGap)
    + logDirichletPdfCounts(pSub,nSub);
}

LogProb MutatorCounts::logLikelihood (const MutatorParams& params) const {
  Assert (nLen.size() == params.pLen.size(), "Length mismatch");
  double ll = 0;
  for (size_t i = 0; i < nLen.size(); ++i)
    ll += nLen[i] * log(params.pLen[i]);
  return nDelOpen * log(params.pDelOpen)
    + nTanDup * log(params.pTanDup)
    + nNoGap * log(params.pNoGap())
    + nTransition() * log(params.pTransition)
    + nTransversion() * log(params.pTransversion)
    + nMatch() * log(params.pMatch())
    + ll;
}
