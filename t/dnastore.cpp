#include <cstdlib>
#include <stdexcept>
#include <boost/program_options.hpp>

#include "../src/vguard.h"
#include "../src/logger.h"

using namespace std;

namespace po = boost::program_options;

typedef unsigned long long Kmer;
typedef unsigned char Base;
typedef int Pos;

const char* alphabet = "ACGT";
inline char baseToChar (Base base) {
  return alphabet[base & 3];
}

inline Base getBase (Kmer kmer, Pos pos) {
  return (kmer >> ((pos - 1) << 1)) & 3;
}

inline string kmerToString (Kmer kmer, Pos len) {
  string s (len, '*');
  for (Pos i = 1; i <= len; ++i)
    s[i-1] = baseToChar(getBase(kmer,i));
  return s;
}

inline bool isTransition (Base x, Base y) {
  return (x & 1) == (y & 1);
}

inline bool isComplement (Base x, Base y) {
  return (x ^ 3) == y;
}

inline Kmer kmerMask (Pos len) {
  // 4^len - 1 = 2^(2*len) - 1 = (1 << (2*len)) - 1 = (1 << (len << 1)) - 1
  return (((Kmer) 1) << (len << 1)) - 1;
}

inline Kmer kmerSubstring (Kmer kmer, Pos start, Pos len) {
  return (kmer >> ((start - 1) << 1)) & kmerMask (len);
}

inline Kmer kmerRevComp (Kmer kmer, Pos len) {
  Kmer rc = 0;
  for (Pos i = 1; i <= len; ++i)
    rc = (rc << 2) | getBase(kmer,i);
  return rc;
}

struct EditDistanceMatrix {
  typedef short Score;
  const Pos len, minMatchLen;
  const Score editDistanceThreshold;
  vguard<vguard<Score> > editScore;
  EditDistanceMatrix (Pos len, Pos minMatchLen, Score editDistanceThreshold)
    : len (len),
      minMatchLen (minMatchLen),
      editDistanceThreshold (editDistanceThreshold),
      editScore (len+1, vguard<Score> (len+1, 0))
  { }
  static inline Score hammingEditScore (Base x, Base y) {
    return x == y ? 0 : 1;
  }
  inline bool pastThreshold (Kmer seq) {
    if (minMatchLen > len)
      return false;
    for (Pos i = len - minMatchLen + 1; i <= len; ++i)
      for (Pos j = 1; j < i; ++j) {
	Score sc = min (editScore[i-1][j-1] + hammingEditScore (getBase(seq,i), getBase(seq,j)),
			editScore[i-1][j] + 1);
	if (j > 1)
	  sc = min (sc, (Score) (editScore[i][j-1] + 1));
	if (i == len) {
	  if (sc < editDistanceThreshold)
	    return true;
	} else
	  editScore[i][j] = sc;
      }
    return false;
  }
};

struct FoldMatrix {
  typedef unsigned short Score;
  const Pos len;
  const Score threshold;
  Pos minBetween;
  vguard<vguard<Score> > basepairs;
  FoldMatrix (Pos len, Score threshold)
    : len (len),
      threshold (threshold),
      basepairs (len+1, vguard<Score> (len+1, 0)),
      minBetween (2)
  { }
  inline Score basepairScore (Base x, Base y) const {
    return isComplement(x,y) ? 1 : 0;
  }
  inline bool pastThreshold (Kmer seq) {
    for (Pos i = len - 1 - minBetween; i >= 1; --i)
      for (Pos j = i + minBetween + 1; j <= len; ++j) {
	Score sc = max ((Score) (basepairScore (getBase(seq,i), getBase(seq,j)) + basepairs[i+1][j-1]),
			max (basepairs[i+1][j], basepairs[i][j-1]));
	for (Pos k = i + minBetween + 2; k <= j - minBetween - 2; ++k)
	  sc = max (sc, (Score) (basepairs[i][k] + basepairs[k+1][j]));
	if (sc >= threshold)
	  return true;
	basepairs[i][j] = sc;
      }
    return false;
  }
};

inline bool hasExactTandemRepeat (Kmer seq, Pos len, Pos maxRepeatLen) {
  for (Pos repeatLen = 1; repeatLen <= maxRepeatLen; ++repeatLen)
    for (Pos i = 1; i <= len - repeatLen; ++i)
      if (kmerSubstring (seq, i, repeatLen) == kmerSubstring (seq, i + repeatLen, repeatLen))
	return true;
  return false;
}

inline bool tooManyMismatches (Kmer seq1, Kmer seq2, Pos len, unsigned int maxMismatches) {
  unsigned int mismatches = 0;
  for (Pos i = 1; i < len; ++i)
    if (getBase(seq1,i) != getBase(seq2,i))
      if (++mismatches > maxMismatches)
	return true;
  return false;
}

inline bool hasApproxTandemRepeat (Kmer seq, Pos len, Pos repeatLen, double minIdentity) {
  if (repeatLen < 1)
    return false;
  const unsigned int maxMismatches = (unsigned int) ((1. - minIdentity) * repeatLen);
  for (Pos i = 1; i <= len - repeatLen; ++i)
    if (!tooManyMismatches (kmerSubstring (seq, i, repeatLen),
			    kmerSubstring (seq, i + repeatLen, repeatLen),
			    repeatLen,
			    maxMismatches))
      return true;
  return false;
}

inline bool hasExactInvertedRepeat (Kmer seq, Pos len, Pos repeatLen) {
  const Kmer rc = kmerRevComp (seq, len);
  for (Pos i = len - repeatLen*2; i > 0; --i) {
    const Kmer invRep = kmerSubstring (rc, len - i + 1, repeatLen);
    const Pos jMin = i + repeatLen;
    for (Pos j = len - repeatLen + 1; j >= jMin; --j)
      if (invRep == kmerSubstring (seq, j, repeatLen))
	return true;
  }
  return false;
}

int main (int argc, char** argv) {

  try {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "display this help message")
      ("kmerlen,k", po::value<int>()->default_value(15), "length of k-mers in de Bruijn graph")
      ("invrep,i", po::value<int>()->default_value(5), "minimum length at which to reject exact nonlocal inverted repeats")
      ("tandem,t", po::value<int>()->default_value(4), "maximum length of exact local tandem repeats")
      ("approx,a", po::value<double>()->default_value(.5), "minimum identity level at which to reject approximate local tandem repeats")
      ("remotelen,l", po::value<int>()->default_value(5), "minimum length at which to reject approximate nonlocal tandem repeats")
      ("remotedist,r", po::value<int>()->default_value(3), "maximum Levenshtein edit distance at which to reject approximate nonlocal tandem repeats")
      ("basepairs,b", po::value<int>()->default_value(6), "minimum number of nested complementary base-pairs at which to reject k-mer")
      ("debug,d", "print all sequences and tests")
      ("verbose,v", po::value<int>()->default_value(1), "verbosity level")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    logger.setVerbose (vm.at("verbose").as<int>());
    const bool debug = vm.count("debug");
    
    if (vm.count("help")) {
      cout << desc << "\n";
      return 1;
    }

    const Pos len = vm.at("kmerlen").as<int>();

    const Pos maxExactTandemRepeatLen = vm.at("tandem").as<int>();
    const Pos minApproxTandemRepeatLen = maxExactTandemRepeatLen;
    const double minApproxTandemRepeatIdentity = vm.at("approx").as<double>();

    const Pos minInvRepLen = vm.at("invrep").as<int>();

    EditDistanceMatrix edit (len, vm.at("remotelen").as<int>(), vm.at("remotedist").as<int>());
    FoldMatrix fold (len, vm.at("basepairs").as<int>());

    const Kmer maxKmer = kmerMask(len);

    ProgressLog (plog, 1);
    plog.initProgress ("Iterating over %d-mers", len);

    for (Kmer kmer = 0; kmer <= maxKmer; ++kmer) {
      plog.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
      if (debug) {

	cout << kmerToString(kmer,len)
	     << " exactLocalTandem=" << hasExactTandemRepeat(kmer,len,maxExactTandemRepeatLen)
	     << " approxLocalTandem=" << hasApproxTandemRepeat(kmer,len,minApproxTandemRepeatLen,minApproxTandemRepeatIdentity)
	     << " exactRemoteInvRep=" << hasExactInvertedRepeat(kmer,len,minInvRepLen)
	     << " approxRemoteTandem=" << edit.pastThreshold(kmer)
	     << " tooManyBasepairs=" << fold.pastThreshold(kmer)
	     << endl;

      } else {

	plog.logProgress (kmer / (double) maxKmer, "sequence %llu/%llu", kmer, maxKmer);
	
	if (!hasExactTandemRepeat(kmer,len,maxExactTandemRepeatLen)
	    && !hasApproxTandemRepeat(kmer,len,minApproxTandemRepeatLen,minApproxTandemRepeatIdentity)
	    && !hasExactInvertedRepeat(kmer,len,minInvRepLen)
	    && !edit.pastThreshold(kmer)
	    && !fold.pastThreshold(kmer))
	  cout << kmerToString(kmer,len) << endl;
      }
    }

  } catch (const std::exception& e) {
    cerr << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}
