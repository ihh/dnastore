#include <fstream>
#include "mutator.h"
#include "jsonutil.h"

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
