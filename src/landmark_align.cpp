// SPDX-License-Identifier: GPL-3.0-or-later
// Fast landmark-based alignment kernels for MALDI-TOF spectra
// Provides:
//  - cpp_match_anchors: greedy 1:1 matching between spectrum peaks and reference anchors
//  - cpp_warp_mz_piecewise: piecewise-linear warping of mass axis through matched anchors
//
// These functions are called from R via Rcpp. See R/alignSpectra_landmark_cpp.R
// and R/SignalProcessingUltra.R for usage.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

/*
  Match sorted spectrum peaks to sorted reference anchors within a symmetric window.

  Arguments
    peaks: numeric vector of spectrum peak m/z (ascending).
    ref  : numeric vector of reference anchor m/z (ascending).
    tol  : tolerance value; interpreted as Da when ppm == false; as ppm when ppm == true.
    ppm  : if true, tolerance is in ppm (per reference mass).

  Returns
    List with:
      src: matched spectrum peak m/z values (source)
      dst: corresponding reference anchor m/z values (destination)
*/
// [[Rcpp::export]]
List cpp_match_anchors(NumericVector peaks,
                       NumericVector ref,
                       double tol,
                       bool ppm) {
  const int n = peaks.size();
  const int m = ref.size();

  std::vector<double> src;
  std::vector<double> dst;
  src.reserve(std::min(n, m));
  dst.reserve(std::min(n, m));

  int i = 0, j = 0;

  // Greedy merge scan with nearest-anchor tie-break
  while (i < n && j < m) {
    const double g = ref[j];
    const double t = ppm ? g * tol * 1e-6 : tol; // convert ppm to Da window
    const double lo = g - t;
    const double hi = g + t;
    const double p = peaks[i];

    if (p < lo) {
      ++i; // spectrum peak is left of current window
      continue;
    }
    if (p > hi) {
      ++j; // window lies left of spectrum peak
      continue;
    }

    // p within [lo, hi]; if next ref is also within, choose the nearer one
    if (j + 1 < m) {
      const double g2 = ref[j + 1];
      const double t2 = ppm ? g2 * tol * 1e-6 : tol;
      const double lo2 = g2 - t2;
      const double hi2 = g2 + t2;
      if (p >= lo2 && p <= hi2) {
        const double d1 = std::abs(p - g);
        const double d2 = std::abs(p - g2);
        if (d2 < d1) {
          ++j; // consider the closer anchor
          continue;
        }
      }
    }

    // Record 1:1 match and advance both pointers
    src.push_back(p);
    dst.push_back(g);
    ++i;
    ++j;
  }

  return List::create(_["src"] = wrap(src),
                      _["dst"] = wrap(dst));
}

/*
  Warp a mass axis by piecewise-linear mapping defined by matched anchor pairs.

  Arguments
    mz         : numeric vector of original m/z values (arbitrary order, typically ascending).
    src_anchor : matched source anchor m/z (ascending, from spectrum).
    dst_anchor : matched destination anchor m/z (ascending, from reference).

  Returns
    Numeric vector with warped m/z values. If:
      - no anchors: returns input unchanged,
      - one anchor: applies constant shift to all m/z,
      - two or more anchors: applies y = a_j * x + b_j within each segment and
        extrapolates beyond the outer anchors using the end-segment slopes.
*/
// [[Rcpp::export]]
NumericVector cpp_warp_mz_piecewise(NumericVector mz,
                                    NumericVector src_anchor,
                                    NumericVector dst_anchor) {
  const int n = mz.size();
  const int k = src_anchor.size();

  NumericVector out = clone(mz); // default: identity

  if (k == 0) {
    return out;
  }

  if (k == 1) {
    const double shift = dst_anchor[0] - src_anchor[0];
    for (int i = 0; i < n; ++i) {
      out[i] = mz[i] + shift;
    }
    return out;
  }

  // Precompute segment coefficients for y = a * x + b on each [src[j], src[j+1]]
  std::vector<double> a(k - 1), b(k - 1);
  for (int j = 0; j < k - 1; ++j) {
    const double x1 = src_anchor[j];
    const double x2 = src_anchor[j + 1];
    const double y1 = dst_anchor[j];
    const double y2 = dst_anchor[j + 1];
    const double denom = (x2 - x1);

    if (std::abs(denom) < 1e-12) {
      // Degenerate segment: fall back to unit slope and shift
      a[j] = 1.0;
      b[j] = y1 - x1;
    } else {
      a[j] = (y2 - y1) / denom;
      b[j] = y1 - a[j] * x1;
    }
  }

  const double aL = a.front();
  const double bL = b.front();
  const double aR = a.back();
  const double bR = b.back();

  // Walk along mz and apply the appropriate segment (two-pointer for segments)
  int seg = 0;
  for (int i = 0; i < n; ++i) {
    const double x = mz[i];

    if (x <= src_anchor[0]) {
      out[i] = aL * x + bL;     // left extrapolation
      continue;
    }
    if (x >= src_anchor[k - 1]) {
      out[i] = aR * x + bR;     // right extrapolation
      continue;
    }
    // advance seg so x in [src_anchor[seg], src_anchor[seg+1]]
    while (seg + 1 < k && x > src_anchor[seg + 1]) {
      ++seg;
    }
    out[i] = a[seg] * x + b[seg];
  }

  return out;
}
