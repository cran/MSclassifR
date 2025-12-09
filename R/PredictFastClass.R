#' Fast class prediction from peak lists using linear regressions
#'
#' Builds a sample-by-m/z matrix from a list of MALDIquant MassPeaks and predicts
#' the class of each spectrum by fitting, for each class, a linear regression of
#' the spectrum’s intensities on the training spectra of that class. The class
#' minimizing the AIC is selected as the predicted label. In parallel, an F-test
#' p-value is computed per class to quantify how unlikely the spectrum is to
#' belong to the training database; the minimum across classes is returned as
#' `p_not_in_DB`. The peak-to-m/z matching is done in C++ via
#' [build_X_from_peaks_fast()] for speed.
#'
#' @param peaks List of MALDIquant::MassPeaks objects to classify (one per spectrum).
#'   Each element must expose `@mass` (numeric m/z) and `@intensity` (numeric) of
#'   the same length. Names/metaData are used to populate the `name` column.
#' @param mod_peaks Numeric training matrix of dimension n_train x p (rows =
#'   spectra, columns = m/z features) used as regressors per class. Column names
#'   must be m/z values (character) and must include all m/z requested in `moz`.
#' @param Y_mod_peaks Factor of length n_train giving the class labels for rows of
#'   `mod_peaks`.
#' @param moz Either "ALL" or a numeric vector of target m/z. If "ALL" (default),
#'   the column names of `mod_peaks` are used. Otherwise, the provided m/z are used
#'   (they must all be present among the column names of `mod_peaks`).
#' @param tolerance Numeric (Da). A target m/z is matched to the nearest peak only
#'   if the absolute difference is <= `tolerance`. Default 6.
#' @param normalizeFun Logical; if TRUE, per-spectrum max normalization is applied
#'   after matching (i.e., each row of the new matrix is divided by its maximum).
#'   Default TRUE.
#' @param noMatch Numeric; intensity value inserted when no peak is matched for a
#'   given target m/z. Default 0.
#' @param chunk_size Integer; rows per block when building the new matrix from
#'   `peaks` (passed to [build_X_from_peaks_fast()], if used). Default 2000.
#' @param ncores Integer; number of cores to use when building the new matrix from
#'   `peaks` (R side). Default 1.
#' @param verbose Logical; print progress messages. Default FALSE.
#'
#' @return A data.frame with columns:
#'   - name: spectrum name (from MassPeaks metaData fullName/file if available).
#'   - p_not_in_DB: minimum F-test p-value across classes (smaller suggests the
#'     spectrum matches the training database; larger suggests “not in DB”).
#'   - pred_cat: predicted class (label with smallest AIC).
#'
#' @details
#' - Matrix building: [build_X_from_peaks_fast()] maps each spectrum in `peaks`
#'   to the target m/z grid with nearest-within-tolerance matching (C++). If
#'   `normalizeFun = TRUE`, each row is divided by its maximum (guarded to avoid
#'   division by zero). Spectra with initially no matches are retried with a
#'   slightly increased tolerance (internal bumping).
#' - Alignment to training: columns of the new matrix must align to `mod_peaks`.
#'   The function stops if any requested m/z is missing from `mod_peaks`.
#' - Per-class regression: for each class k, it regresses the new spectrum’s
#'   intensities on the columns of `mod_peaks` belonging to class k (after
#'   removing entries where the new spectrum is non-finite). If the number of
#'   training spectra exceeds the number of non-missing points in the spectrum,
#'   a random subset of columns (size = length(non-missing) - 1) is used to avoid
#'   singular fits. Fitting is done via stats::lm.fit for speed.
#' - Selection and scores: `pred_cat` is the class with smallest AIC across fitted
#'   models. For each class, an F-test p-value is computed from the model summary;
#'   `p_not_in_DB` is the minimum across classes (1 if a class model fails).
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("MALDIquant", quietly = TRUE)) {
#'   set.seed(1)
#'   # Create a small training set (mod_peaks) with 2 classes
#'   p <- 6
#'   moz <- as.character(round(seq(1000, 1500, length.out = p), 2))
#'   mod_peaks <- rbind(
#'     matrix(runif(5 * p, 0, 1), nrow = 5, dimnames = list(NULL, moz)),
#'     matrix(runif(5 * p, 0, 1), nrow = 5, dimnames = list(NULL, moz))
#'   )
#'   Y_mod <- factor(rep(c("A", "B"), each = 5))
#'
#'   # Two spectra to classify: generate MassPeaks near moz
#'   mk_peaks <- function(shift = 0) {
#'     MALDIquant::createMassPeaks(
#'       mass = as.numeric(moz) + rnorm(length(moz), shift, 0.2),
#'       intensity = runif(length(moz), 10, 100)
#'     )
#'   }
#'   peaks <- list(mk_peaks(0.1), mk_peaks(-0.1))
#'
#'   res <- PredictFastClass(
#'     peaks = peaks,
#'     mod_peaks = mod_peaks,
#'     Y_mod_peaks = Y_mod,
#'     moz = "ALL",
#'     tolerance = 1,
#'     normalizeFun = TRUE
#'   )
#'   res
#' }
#' }
#'
#' @seealso build_X_from_peaks_fast; MALDIquant::createMassPeaks; stats::lm.fit
#' @export
PredictFastClass <- function(peaks,
                                  mod_peaks,
                                  Y_mod_peaks,
                                  moz = "ALL",
                                  tolerance = 6,
                                  normalizeFun = TRUE,
                                  noMatch = 0,
                                  chunk_size = 2000L,
                                  ncores = 1L,
                                  verbose = FALSE) {
  Y <- factor(Y_mod_peaks)
  if (identical(moz, "ALL")) moz <- colnames(mod_peaks)
  moz_num <- sort(unique(as.numeric(moz)))

  # Build new X once
  Xnew <- build_X_from_peaks_fast(peaks, moz_num, tolerance, normalize = normalizeFun,
                                  noMatch = noMatch, bump_if_empty = TRUE, toleranceStep = 2)

  # Align columns of Xnew to mod_peaks
  idx_cols <- match(colnames(Xnew), colnames(mod_peaks))
  if (anyNA(idx_cols)) stop("mod_peaks is missing some requested m/z columns")
  IntM <- t(mod_peaks[, idx_cols, drop = FALSE])  # p' x n_train
  desig <- stats::model.matrix(~ Y - 1)

  n <- nrow(Xnew)
  name <- vapply(seq_along(peaks), function(i) {
    pk <- peaks[[i]]
    nm <- pk@metaData$fullName
    if (is.null(nm) || length(nm) == 0L || !nzchar(as.character(nm)[1])) {
      nm <- pk@metaData$file
    }
    if (is.null(nm) || length(nm) == 0L || !nzchar(as.character(nm)[1])) {
      # fallback to list name or generated id
      nm <- if (!is.null(names(peaks)) && length(names(peaks)) >= i &&
                nzchar(names(peaks)[i])) names(peaks)[i] else paste0("spec_", i)
    }
    as.character(nm)[1]
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
  pvalAIC <- matrix(NA_real_, nrow = n, ncol = nlevels(Y), dimnames = list(NULL, levels(Y)))
  pvalF <- matrix(NA_real_, nrow = n, ncol = nlevels(Y), dimnames = list(NULL, levels(Y)))

  for (i in seq_len(n)) {
    newXIC <- Xnew[i, ]
    ok <- is.finite(newXIC)
    if (sum(ok) < 2L) { pvalF[i, ] <- 1; next }
    yy <- (newXIC[ok]) * 1e10

    for (k in seq_len(nlevels(Y))) {
      cols <- which(desig[, k] == 1)
      if (length(cols) == 0L) { pvalF[i, k] <- 1; next }
      Bb <- IntM[ok, cols, drop = FALSE]
      Bb[!is.finite(Bb)] <- 0
      if (ncol(Bb) >= length(yy)) {
        set.seed(1L)
        Bb <- Bb[, sample.int(ncol(Bb), max(1L, length(yy) - 1L)), drop = FALSE]
      }
      fit <- try(stats::lm.fit(x = cbind(1, Bb), y = yy), silent = TRUE)
      if (inherits(fit, "try-error") || fit$df.residual <= 0L) {
        pvalF[i, k] <- 1
        next
      }
      rss <- sum(fit$residuals^2)
      kpar <- ncol(Bb) + 1L
      pvalAIC[i, k] <- length(yy) * log(rss / length(yy)) + 2 * kpar
      ss_tot <- sum((yy - mean(yy))^2)
      ss_reg <- ss_tot - rss
      df1 <- kpar - 1L; df2 <- fit$df.residual
      if (df1 > 0 && df2 > 0 && rss > 0) {
        Fstat <- (ss_reg / df1) / (rss / df2)
        pvalF[i, k] <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)
      } else {
        pvalF[i, k] <- 1
      }
    }
  }

  pred_idx <- apply(pvalAIC, 1L, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
  pred_cat <- levels(Y)[pred_idx]
  min_p <- apply(pvalF, 1L, function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else 1 })
  data.frame(name = name, p_not_in_DB = min_p, pred_cat = pred_cat, check.names = FALSE)
}
