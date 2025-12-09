// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Compute Kruskal–Wallis p-values column-wise.
// grp: integer group indices 0..(G-1), length n. NA in X skipped.
// [[Rcpp::export]]
NumericVector kruskal_cols_cpp(NumericMatrix X, IntegerVector grp) {
  const int n = X.nrow();
  const int p = X.ncol();
  int G = 0;
  for (int i = 0; i < grp.size(); ++i) if (grp[i] > G) G = grp[i];
  G += 1;

  NumericVector out(p, NA_REAL);
  if (n == 0 || p == 0 || G < 2) return out;

  for (int j = 0; j < p; ++j) {
    // collect finite values and group ids
    std::vector<double> vals;
    std::vector<int> gids;
    vals.reserve(n); gids.reserve(n);
    for (int i = 0; i < n; ++i) {
      double v = X(i, j);
      if (R_finite(v)) { vals.push_back(v); gids.push_back(grp[i]); }
    }
    const int m = (int)vals.size();
    if (m < 2) { out[j] = NA_REAL; continue; }

    // ranks with average ties
    std::vector<int> ord(m);
    for (int i = 0; i < m; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return vals[a] < vals[b]; });

    std::vector<double> ranks(m);
    double tie_correction_num = 0.0; // sum(t^3 - t)
    int r = 1;
    for (int k = 0; k < m; ) {
      int k2 = k + 1;
      while (k2 < m && vals[ord[k2]] == vals[ord[k]]) ++k2;
      const int t = k2 - k;
      const double avg_rank = (r + (r + t - 1)) / 2.0;
      for (int u = k; u < k2; ++u) ranks[ord[u]] = avg_rank;
      if (t > 1) tie_correction_num += (double)t * (t * t - 1.0);
      r += t;
      k = k2;
    }

    // sum ranks per group
    std::vector<double> Rg(G, 0.0);
    std::vector<int> ng(G, 0);
    for (int i = 0; i < m; ++i) {
      int g = gids[i];
      if (g >= 0 && g < G) { Rg[g] += ranks[i]; ng[g] += 1; }
    }

    // effective groups
    int Gm = 0;
    for (int g = 0; g < G; ++g) if (ng[g] > 0) Gm++;
    if (Gm < 2) { out[j] = NA_REAL; continue; }

    const double N = (double)m;
    double sum_term = 0.0;
    for (int g = 0; g < G; ++g) if (ng[g] > 0) sum_term += (Rg[g] * Rg[g]) / (double)ng[g];

    double H = (12.0 / (N * (N + 1.0))) * sum_term - 3.0 * (N + 1.0);

    // tie correction
    double T = 1.0;
    if (tie_correction_num > 0.0) {
      T = 1.0 - tie_correction_num / (N * N * N - N);
      if (T <= 0.0) T = 1.0;
    }
    H /= T;

    // upper-tail p-value, df = Gm - 1
    out[j] = R::pchisq(H, (double)(Gm - 1), /*lower_tail*/ 0, /*log_p*/ 0);
  }
  return out;
}

// Compute one-way ANOVA p-values column-wise.
// grp: integer group indices 0..(G-1). NA in X skipped.
// [[Rcpp::export]]
NumericVector anova_cols_cpp(NumericMatrix X, IntegerVector grp) {
  const int n = X.nrow();
  const int p = X.ncol();
  int G = 0;
  for (int i = 0; i < grp.size(); ++i) if (grp[i] > G) G = grp[i];
  G += 1;

  NumericVector out(p, NA_REAL);
  if (n == 0 || p == 0 || G < 2) return out;

  for (int j = 0; j < p; ++j) {
    std::vector<double> sumg(G, 0.0), sumsqg(G, 0.0);
    std::vector<int> ng(G, 0);
    double sum_all = 0.0;
    int n_all = 0;

    for (int i = 0; i < n; ++i) {
      double v = X(i, j);
      if (R_finite(v)) {
        int g = grp[i];
        sumg[g]   += v;
        sumsqg[g] += v * v;
        ng[g]     += 1;
        sum_all   += v;
        n_all     += 1;
      }
    }
    if (n_all < 2) { out[j] = NA_REAL; continue; }

    double mean_all = sum_all / (double)n_all;

    double SSB = 0.0, SSE = 0.0;
    int Gm = 0;
    for (int g = 0; g < G; ++g) if (ng[g] > 0) {
      double mean_g = sumg[g] / (double)ng[g];
      SSB += (double)ng[g] * (mean_g - mean_all) * (mean_g - mean_all);
      SSE += sumsqg[g] - (sumg[g] * sumg[g]) / (double)ng[g];
      Gm++;
    }
    if (Gm < 2 || (n_all - Gm) <= 0) { out[j] = NA_REAL; continue; }

    double df1 = (double)(Gm - 1);
    double df2 = (double)(n_all - Gm);
    double MSB = SSB / df1;
    double MSE = SSE / df2;
    if (MSE <= 0.0) { out[j] = 1.0; continue; }

    double F = MSB / MSE;
    out[j] = R::pf(F, df1, df2, /*lower_tail*/ 0, /*log_p*/ 0);
  }
  return out;
}
