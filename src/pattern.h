#ifndef PATTERN_INCLUDED
#define PATTERN_INCLUDED

#include <set>
#include "kmer.h"
#include "logger.h"

inline bool endsWithMotif (Kmer seq, Pos len, const set<KmerLen>& motif, const char* desc = NULL) {
  for (const auto& kl: motif)
    if (kmerSub(seq,1,kl.len) == kl.kmer) {
      if (desc)
	LogThisAt(4,"Rejecting " << kmerString(seq,len) << " because it ends with " << kmerString(kl) << " (" << desc << ")" << endl);
      return true;
    }
  return false;
}

inline bool hasExactTandemRepeat (Kmer seq, Pos len, Pos maxRepeatLen) {
  for (Pos repeatLen = 1; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = len - 2*repeatLen + 1; i >= 1; --i)
      if (kmerSub (seq, i, repeatLen) == kmerSub (seq, i + repeatLen, repeatLen)) {
	const int logLevel = max(5,8-repeatLen);
	LogThisAt(logLevel,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,i+repeatLen,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (" << (repeatLen == 1 ? "repeated base" : "exact tandem repeat") << ")" << endl);
	return true;
      }
  return false;
}

inline bool hasExactLocalInvertedRepeat (Kmer seq, Pos len, Pos minRepeatLen, Pos maxRepeatLen) {
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos repeatLen = minRepeatLen; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = len - 2*repeatLen + 1; i >= 1; --i) {
      const Kmer invRep = kmerSub (rc, len - i + 1, repeatLen);
      if (invRep == kmerSub (seq, i + repeatLen, repeatLen)) {
	const int logLevel = max(5,8-repeatLen);
	LogThisAt(logLevel,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,i+repeatLen,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (palindrome)" << endl);
	return true;
      }
    }
  return false;
}

inline bool hasExactNonlocalInvertedRepeat (Kmer seq, Pos len, Pos repeatLen, Pos minSeparation) {
  if (repeatLen <= 0)
    return false;
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos i = len - repeatLen*2 - minSeparation; i > 0; --i) {
    const Kmer invRep = kmerSub (rc, len - i + 1, repeatLen);
    const Pos jMin = i + repeatLen + minSeparation;
    for (Pos j = len - repeatLen + 1; j >= jMin; --j)
      if (invRep == kmerSub (seq, j, repeatLen)) {
	LogThisAt(4,"Rejecting " << kmerString(seq,len) << " because " << kmerSubAt(seq,j,repeatLen,len) << " matches " << kmerSubAt(seq,i,repeatLen,len) << " (exact inverted repeat)" << endl);
	return true;
      }
  }
  return false;
}

#endif /* PATTERN_INCLUDED */
