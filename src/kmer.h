#ifndef KMER_INCLUDED
#define KMER_INCLUDED

#include <cmath>
#include <string>
#include "vguard.h"
#include "util.h"

using namespace std;

typedef unsigned long long Kmer;
typedef unsigned short Base;
typedef int Pos;

struct KmerLen {
  Kmer kmer;
  Pos len;
  KmerLen (Kmer kmer, Pos len) : kmer(kmer), len(len) { }
  bool operator< (const KmerLen& kl) const {
    return len < kl.len || (len == kl.len && kmer < kl.kmer);
  }
};

extern const char* dnaAlphabet;  // ACGT
extern string dnaAlphabetString;

#define AdenineBase  0
#define CytosineBase 1
#define GuanineBase  2
#define ThymineBase  3

#define BaseMask       3
#define PyrimidineMask 1
#define CarbonylMask   2  /* G,T have one more carbonyl side-group than (respectively) A,C */

inline char baseToChar (Base base) {
  return dnaAlphabet[base & BaseMask];
}

inline Base charToBase (char c) {
  const char* s = strchr (dnaAlphabet, toupper(c));
  Require (s != NULL, "%c is not a nucleotide character", c);
  return (Base) (s - dnaAlphabet);
}

inline Base getBase (Kmer kmer, Pos pos) {
  return (kmer >> ((pos - 1) << 1)) & BaseMask;
}

inline Kmer setBase (Kmer kmer, Pos pos, Base base) {
  const int shift = (pos - 1) << 1;
  return (kmer & (((Kmer) -1) ^ (BaseMask << shift))) | (((Kmer) base) << shift);
}

inline Base complementBase (Base b) {
  return BaseMask - b;
}

inline Kmer makeTransition (Kmer kmer, Pos pos) {
  return kmer ^ (CarbonylMask << ((pos - 1) << 1));
}

inline string kmerString (Kmer kmer, Pos len) {
  string s (len, '*');
  for (Pos i = 1; i <= len; ++i)
    s[len-i] = baseToChar(getBase(kmer,i));
  return s;
}

inline string kmerString (KmerLen kl) {
  return kmerString (kl.kmer, kl.len);
}

inline Kmer stringToKmer (const string& s) {
  Kmer kmer = 0;
  for (Pos i = 0; i <= s.size(); ++i)
    kmer = setBase (kmer, s.size() - i, charToBase (s[i]));
  return kmer;
}

inline Kmer stringToKmer (const char* s) {
  return stringToKmer (string (s));
}

inline bool isTransition (Base x, Base y) {
  return x != y && (x & PyrimidineMask) == (y & PyrimidineMask);
}

inline bool isTransversion (Base x, Base y) {
  return x != y && (x & PyrimidineMask) != (y & PyrimidineMask);
}

inline bool isComplement (Base x, Base y) {
  return y == complementBase(x);
}

inline bool isGC (Base x) {
  return x == GuanineBase || x == CytosineBase;
}

inline double gcContent (Kmer kmer, Pos len) {
  int gc = 0;
  for (Pos pos = 1; pos <= len; ++pos)
    if (isGC (getBase (kmer, pos)))
      ++gc;
  return gc / (double) len;
}

inline double gcNonuniformity (Kmer kmer, Pos len) {
  const double gc = gcContent (kmer, len);
  return abs (gc - .5);
}

inline double kmerEntropy (Kmer kmer, Pos len) {
  vguard<int> freq (4);
  for (Pos pos = 1; pos <= len; ++pos)
    ++freq[getBase(kmer,pos)];
  double S = 0;
  for (Base b = 0; b < 3; ++b)
    if (freq[b])
      S -= freq[b] * log (freq[b] / (double) len);
  return S / log(2);
}

inline size_t kmerHammingDistance (Kmer a, Kmer b, Pos len) {
  size_t d = 0;
  for (Pos pos = 1; pos <= len; ++pos)
    if (getBase(a,pos) != getBase(b,pos))
      ++d;
  return d;
}

inline Kmer kmerMask (Pos len) {
  // 4^len - 1 = 2^(2*len) - 1 = (1 << (2*len)) - 1 = (1 << (len << 1)) - 1
  return (((Kmer) 1) << (len << 1)) - 1;
}

inline Kmer kmerSub (Kmer kmer, Pos start, Pos len) {
  return (kmer >> ((start - 1) << 1)) & kmerMask (len);
}

inline string kmerSubstring (Kmer kmer, Pos start, Pos len) {
  return kmerString (kmerSub (kmer, start, len), len);
}

inline int kmerLeftCoord (Pos pos, Pos len) {
  return len - pos + 1;
}

inline string kmerSubCoords (Pos start, Pos len, Pos kmerLen) {
  return string("[") + to_string(kmerLeftCoord(start+len-1,kmerLen)) + (len > 1 ? (string("..") + to_string(kmerLeftCoord(start,kmerLen))) : string()) + "]";
}

inline string kmerSubAt (Kmer kmer, Pos start, Pos len, Pos kmerLen) {
  return kmerString (kmerSub (kmer, start, len), len) + kmerSubCoords(start,len,kmerLen);
}

inline Kmer kmerRevComp (Kmer kmer, Pos len) {
  Kmer rc = 0;
  for (Pos i = 1; i <= len; ++i)
    rc = (rc << 2) | complementBase (getBase(kmer,i));
  return rc;
}

struct EdgeVector : vguard<Kmer> {
  EdgeVector() : vguard<Kmer>(4) { }
};

#endif /* KMER_INCLUDED */
