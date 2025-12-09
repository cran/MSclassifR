#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Map a single spectrum to moz (nearest within tol). moz must be sorted.
// [[Rcpp::export]]
NumericVector map_spectrum_to_moz_cpp(NumericVector moz,
                                      NumericVector peak_mz,
                                      NumericVector peak_int,
                                      double tol,
                                      double noMatch = 0.0) {
  const int p = moz.size();
  const int k = peak_mz.size();
  NumericVector out(p, noMatch);
  if (p == 0 || k == 0) return out;

  // order peaks by m/z
  std::vector<int> ord(k);
  for (int i = 0; i < k; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return peak_mz[a] < peak_mz[b]; });

  std::vector<double> y(k), v(k);
  y.reserve(k); v.reserve(k);
  for (int i = 0; i < k; ++i) {
    double mzv = peak_mz[ord[i]];
    double iv  = peak_int[ord[i]];
    y[i] = mzv;
    v[i] = iv;
  }

  for (int i = 0; i < p; ++i) {
    double x = moz[i];

    // lower_bound: first y[j] >= x
    auto it = std::lower_bound(y.begin(), y.end(), x);
    int j2 = static_cast<int>(it - y.begin());
    int j1 = j2 - 1;

    double bestDist = std::numeric_limits<double>::infinity();
    double bestVal  = noMatch;

    if (j1 >= 0) {
      double d = std::fabs(x - y[j1]);
      if (d < bestDist) { bestDist = d; bestVal = v[j1]; }
    }
    if (j2 < k) {
      double d = std::fabs(x - y[j2]);
      if (d < bestDist) { bestDist = d; bestVal = v[j2]; }
    }
    if (bestDist <= tol) out[i] = bestVal; // else remains noMatch
  }
  return out;
}

// Build X (n x p) for many spectra. moz must be sorted.
// [[Rcpp::export]]
NumericMatrix build_X_from_peaks_cpp(NumericVector moz,
                                     List mass_list,
                                     List int_list,
                                     double tol,
                                     double noMatch = 0.0,
                                     bool normalize = true) {
  const int n = mass_list.size();
  const int p = moz.size();
  NumericMatrix X(n, p);

  for (int i = 0; i < n; ++i) {
    NumericVector mz_i = mass_list[i];
    NumericVector in_i = int_list[i];

    NumericVector row = map_spectrum_to_moz_cpp(moz, mz_i, in_i, tol, noMatch);

    if (normalize) {
      double rmax = 0.0;
      for (int j = 0; j < p; ++j) {
        double val = row[j];
        if (std::isfinite(val) && val > rmax) rmax = val;
      }
      double denom = (std::isfinite(rmax) && rmax > 0.0) ? rmax : 1.0;
      for (int j = 0; j < p; ++j) {
        double val = row[j];
        X(i, j) = std::isfinite(val) ? (val / denom) : noMatch;
      }
    } else {
      for (int j = 0; j < p; ++j) {
        double val = row[j];
        X(i, j) = std::isfinite(val) ? val : noMatch;
      }
    }
  }
  return X;
}
