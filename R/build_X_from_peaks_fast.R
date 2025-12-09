#' Build a sample-by-m/z intensity matrix from a list of peaks (fast, C++-backed)
#'
#' Converts a list of MALDIquant MassPeaks into a numeric matrix X (rows = samples,
#' columns = target m/z), by matching each target m/z to the nearest peak within
#' a tolerance. If requested, per-row max normalization is applied. Spectra that
#' initially produce no matches can be retried with an increased tolerance.
#' Internally uses Rcpp for speed (binary search per m/z).
#'
#' @param peaks List of MALDIquant::MassPeaks objects (one per sample). Each element
#'   must provide pk@mass (numeric vector of m/z) and pk@intensity (numeric vector
#'   of intensities) of the same length.
#' @param moz Numeric vector of target m/z values (in Da). Will be sorted and uniqued;
#'   the output matrix columns follow this sorted order.
#' @param tolerance Numeric scalar (Da). A target m/z is matched to the nearest
#'   peak only if the absolute difference is <= tolerance. Default 6.
#' @param normalize Logical; if TRUE, per-row max normalization is applied after
#'   matching (i.e., each sample is divided by its maximum non-NA intensity).
#'   Default TRUE.
#' @param noMatch Numeric scalar; intensity value to insert when no peak is matched
#'   for a given target m/z. Default 0.
#' @param bump_if_empty Logical; if TRUE, any spectrum resulting in an all-noMatch
#'   (or all-zero after normalization) row will be retried by increasing the tolerance
#'   in steps of `toleranceStep`, up to `max_bumps` attempts. Default TRUE.
#' @param toleranceStep Numeric scalar (Da); the increment used when bumping the
#'   tolerance for empty rows. Default 2.
#' @param max_bumps Integer; maximum number of bumps when retrying empty rows.
#'   Default 5.
#'
#' @return A numeric matrix X of dimension n x p:
#'   - n = length(peaks)
#'   - p = length(unique(sort(moz)))
#'   Column names are the sorted `moz` coerced to character. Values are intensities
#'   (possibly normalized) or `noMatch` for unmatched positions.
#'
#' @details
#' - Matching: for each target m/z, the nearest peak is chosen if its distance
#'   is <= `tolerance`; otherwise `noMatch` is used. Ties are resolved by nearest
#'   distance via binary search.
#' - Normalization: when `normalize = TRUE`, each row is divided by its maximum
#'   non-NA intensity (guarded to avoid division by zero).
#' - Empty rows: when `bump_if_empty = TRUE`, rows with all `noMatch` (or all zeros
#'   after normalization) are retried with increased tolerance (by `toleranceStep`)
#'   up to `max_bumps` times.
#' - Performance: implemented with C++ helpers (map_spectrum_to_moz_cpp and
#'   build_X_from_peaks_cpp) for speed on large datasets.
#'
#' @examples
#' # Minimal example with synthetic MassPeaks
#' if (requireNamespace("MALDIquant", quietly = TRUE)) {
#'   set.seed(1)
#'   # Two spectra with slightly jittered peaks around 1000, 1500, 2000 Da
#'   mz1 <- c(999.7, 1500.2, 2000.1); int1 <- c(10, 50, 30)
#'   mz2 <- c(1000.3, 1499.8, 2000.4); int2 <- c(12, 60, 28)
#'   p1 <- MALDIquant::createMassPeaks(mass = mz1, intensity = int1)
#'   p2 <- MALDIquant::createMassPeaks(mass = mz2, intensity = int2)
#'   peaks <- list(p1, p2)
#'
#'   # Target m/z grid (unsorted, will be sorted internally)
#'   moz <- c(2000, 1500, 1000)
#'
#'   X <- build_X_from_peaks_fast(peaks, moz, tolerance = 1, normalize = TRUE)
#'   dim(X)
#'   colnames(X)
#'   X
#' }
#'
#' # Typical usage in a pipeline:
#' # spectra <- SignalProcessing(yourSpectra)
#' # peaks   <- MSclassifR::PeakDetection(x = spectra, averageMassSpec = FALSE)
#' # moz     <- c(1000, 1500, 2000)  # from selection or prior knowledge
#' # X       <- build_X_from_peaks_fast(peaks, moz, tolerance = 6, normalize = TRUE)
#' # Then pass X to SelectionVar/SelectionVarStat_fast/LogReg, etc.
#'
#' @seealso MALDIquant::createMassPeaks; internal C++ helpers
#'   map_spectrum_to_moz_cpp and build_X_from_peaks_cpp (not user-facing).
#' @export

build_X_from_peaks_fast <- function(peaks,
                                    moz,
                                    tolerance = 6,
                                    normalize = TRUE,
                                    noMatch = 0,
                                    bump_if_empty = TRUE,
                                    toleranceStep = 2,
                                    max_bumps = 5L) {
  stopifnot(is.list(peaks), length(peaks) > 0L)
  moz <- sort(unique(as.numeric(moz)))
  p <- length(moz)
  n <- length(peaks)

  mass_list <- lapply(peaks, function(pk) pk@mass)
  int_list  <- lapply(peaks, function(pk) pk@intensity)

  X <- build_X_from_peaks_cpp(moz, mass_list, int_list,
                              tol = tolerance, noMatch = noMatch, normalize = normalize)
  colnames(X) <- as.character(moz)

  if (isTRUE(bump_if_empty) && p > 0) {
    # empty row = all noMatch (or all zeros after normalize)
    is_empty <- function(r) {
      if (all(!is.finite(r))) return(TRUE)
      if (all(r == noMatch)) return(TRUE)
      if (normalize && sum(abs(r)) == 0) return(TRUE)
      FALSE
    }
    empty_idx <- which(apply(X, 1L, is_empty))
    if (length(empty_idx)) {
      for (i in empty_idx) {
        tol_i <- tolerance
        tries <- 0L
        xi <- rep(noMatch, p)
        while (tries < max_bumps) {
          tol_i <- tol_i + toleranceStep
          xi <- map_spectrum_to_moz_cpp(moz, mass_list[[i]], int_list[[i]],
                                        tol = tol_i, noMatch = noMatch)
          if (normalize) {
            m <- suppressWarnings(max(xi[is.finite(xi)], na.rm = TRUE))
            if (!is.finite(m) || m <= 0) m <- 1
            xi <- xi / m
          }
          if (!is_empty(xi)) break
          tries <- tries + 1L
        }
        X[i, ] <- xi
      }
    }
  }
  X
}
