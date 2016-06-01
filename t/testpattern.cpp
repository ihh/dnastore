#include "../src/pattern.h"
#include "../src/util.h"

#define TestOK(EXPR) do { if (!(EXPR)) { cout << "Failed: "  #EXPR "\n"; ok = false; } } while (false)

int main (int argc, char** argv) {
  bool ok = true;

  TestOK (hasExactTandemRepeat(stringToKmer("ACGACG"),6,3));
  TestOK (!hasExactTandemRepeat(stringToKmer("ACGACT"),6,3));

  TestOK (hasExactLocalInvertedRepeat(stringToKmer("ACGCGA"),6,1,4));
  TestOK (hasExactLocalInvertedRepeat(stringToKmer("ACGCGA"),6,2,4));
  TestOK (!hasExactLocalInvertedRepeat(stringToKmer("ACGCGA"),6,3,4));
  TestOK (hasExactLocalInvertedRepeat(stringToKmer("ACGCGT"),6,3,4));
  TestOK (!hasExactLocalInvertedRepeat(stringToKmer("ACGCGT"),6,4,4));

  TestOK (!hasExactNonlocalInvertedRepeat (stringToKmer("ACGCGT"),6,3,2));
  TestOK (!hasExactNonlocalInvertedRepeat (stringToKmer("ACGTCGT"),7,3,2));
  TestOK (hasExactNonlocalInvertedRepeat (stringToKmer("ACGTTCGT"),8,3,2));

  cout << (ok ? "ok: pattern recognition works. Kudos to Cayce Pollard" : "not ok: pattern recognition broken. Email Hubertus.Bigend@BlueAnt.com") << endl;

  return EXIT_SUCCESS;
}
