#include <iostream>
#include "logsumexp.h"
#include "util.h"

LogSumExpLookupTable logSumExpLookupTable = LogSumExpLookupTable();

LogSumExpLookupTable::LogSumExpLookupTable() {
  lookup = new double [LOG_SUM_EXP_LOOKUP_ENTRIES];
  int n;
  double x;
  for (n = 0; n < LOG_SUM_EXP_LOOKUP_ENTRIES; ++n) {
    x = n * LOG_SUM_EXP_LOOKUP_PRECISION;
    lookup[n] = log_sum_exp_unary_slow(x);
  }
}

LogSumExpLookupTable::~LogSumExpLookupTable() {
  delete[] lookup;
}

double log_sum_exp_slow (double a, double b) {
  double min, max, diff, ret;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  ret = max + log_sum_exp_unary_slow (diff);
#if defined(NAN_DEBUG)
  if (std::isnan(ret)) {
    cerr << "NaN error in log_sum_exp" << endl;
    throw;
  }
#endif
  return ret;
}

double log_sum_exp_slow (double a, double b, double c) {
  return log_sum_exp_slow (log_sum_exp_slow (a, b), c);
}

double log_sum_exp_slow (double a, double b, double c, double d) {
  return log_sum_exp_slow (log_sum_exp_slow (log_sum_exp_slow (a, b), c), d);
}

double log_sum_exp_unary_slow (double x) {
  return log (1. + exp(-x));
}

void log_accum_exp_slow (double& a, double b) {
  a = log_sum_exp_slow (a, b);
}

double logBetaPdf (double prob, double alpha, double beta) {
  return lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) + (alpha-1)*log(prob) + (beta-1)*log(1-prob);
}

double logDirichletPdf (const vector<double>& prob, const vector<double>& alpha) {
  Assert (prob.size() == alpha.size(), "Dimensionality of Dirichlet counts vector does not match that of probability parameter vector");
  double ld = lgamma (accumulate (alpha.begin(), alpha.end(), 0.));
  for (size_t n = 0; n < prob.size(); ++n)
    ld += (alpha[n] - 1) * log(prob[n]) - lgamma(alpha[n]);
  return ld;
}

double logBetaPdfCounts (double prob, double yesCount, double noCount) {
  return logBetaPdf (prob, yesCount + 1, noCount + 1);
}

double logDirichletPdfCounts (const vector<double>& prob, const vector<double>& count) {
  vector<double> countPlusOne (count);
  for (auto& c : countPlusOne)
    ++c;
  return logDirichletPdf (prob, countPlusOne);
}
