#' Build design matrix X and response Y from peak intensities
#'
#' Constructs a sample-by-peak design matrix (X) and an outcome vector/factor (Y) from
#' peak-intensity input. Accepts either a numeric matrix/data.frame (rows = samples,
#' columns = peaks) or a list of aligned per-sample peak vectors. Optionally applies
#' per-sample max normalization and can return X as a sparse dgCMatrix for memory efficiency.
#'
#' @param peaks Peak data from which to build the design matrix X. Either:
#'   - a numeric matrix/data.frame of intensities with rows = samples and columns = peaks, or
#'   - a list of per-sample numeric vectors (or small tables) aligned to a common set of peaks
#'     (e.g., same names or column order).
#'   Values are assumed non-negative; NAs are allowed and are ignored when computing per-sample maxima.
#'
#' @param labels Outcome/response labels used to build Y. A vector or factor with one entry per sample
#'   (length must equal nrow(peaks) for matrix/data.frame input, or length(peaks) for list input).
#'
#' @param normalize Normalization to apply to each sample’s peak intensities before constructing X.
#'   One of "max" or "none" (matched via match.arg). "max" scales each sample by its maximum
#'   non-NA intensity; "none" applies no scaling.
#'
#' @param sparse Logical; if TRUE, return X as a sparse Matrix::dgCMatrix. If FALSE, return a base R
#'   dense matrix. Default is FALSE.
#'
#' @param name_cols Logical; if TRUE, set column names of X from the m/z values
#'   returned by MALDIquant::intensityMatrix (formatted as "mz_<mz>"). Default is FALSE.
#'   This only applies when peaks is a list of MALDIquant::MassPeaks (or when attr(X, "mass")
#'   is available); otherwise it is ignored. Enabling this can add noticeable overhead for
#'   very wide matrices.
#'
#' @param name_digits Integer scalar; number of decimal digits to use when formatting
#'   m/z values into column names if name_cols = TRUE. Default is 4. Must be a
#'   non-negative integer. Ignored when name_cols = FALSE or when no m/z vector is
#'   available (i.e., for plain matrix/data.frame or generic list inputs).
#'
#' @return A list with:
#'   - X: numeric matrix or Matrix::dgCMatrix of dimension n_samples x n_peaks
#'   - Y: response vector/factor aligned to rows of X (returned as supplied/coerced by the function)
#'
#' @examples
#' data("CitrobacterRKIspectra", "CitrobacterRKImetadata", package = "MSclassifR")
#'
#' spectra <- SignalProcessing(CitrobacterRKIspectra)
#' peaks <- MSclassifR::PeakDetection(x = spectra, averageMassSpec = FALSE)
#'
#' labels <- CitrobacterRKImetadata$Species   # adjust to your label column
#'
#' xy <- build_XY_from_peaks(peaks, labels, normalize = "max", sparse = TRUE)
#'
build_XY_from_peaks <- function(peaks,
                                labels,
                                normalize = c("max", "none"),
                                sparse = FALSE,
                                name_cols = FALSE,
                                name_digits = 4) {
  normalize <- match.arg(normalize)

  # Coerce input to a dense intensity matrix X in the most direct way
  is_masspeaks_list <- is.list(peaks) && length(peaks) > 0L &&
    all(vapply(peaks, inherits, logical(1), what = "MassPeaks"))

  if (is_masspeaks_list) {
    # list of MALDIquant::MassPeaks (directly after PeakDetection)

    # Precompute per-sample maxima cheaply from MassPeaks if needed
    if (normalize == "max") {
      m <- vapply(peaks,
                  function(pk) {
                    v <- MALDIquant::intensity(pk)
                    if (length(v)) max(v, na.rm = TRUE) else 0
                  },
                  numeric(1))
      m[!is.finite(m) | m <= 0] <- 1
    } else {
      m <- NULL
    }

    X <- MALDIquant::intensityMatrix(peaks)  # dense matrix; attr 'mass' holds m/z
    masses <- attr(X, "mass")

    # Optional column naming (skip for speed by default)
    if (name_cols && !is.null(masses)) {
      # Fast fixed-width formatting
      colnames(X) <- paste0("mz_", sprintf(paste0("%.", name_digits, "f"), masses))
    }

    # Normalize and/or convert to sparse with minimal copying
    if (isTRUE(sparse)) {
      X <- Matrix::Matrix(X, sparse = TRUE)
      if (normalize == "max") {
        X <- Matrix::Diagonal(x = 1 / m) %*% X
      }
      if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
    } else {
      if (normalize == "max") {
        X <- sweep(X, 1L, m, "/")
      }
    }

  } else if (is.matrix(peaks)) {
    X <- peaks
    storage.mode(X) <- "double"
    if (normalize == "max") {
      # Fast row maxima for dense matrices
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        m <- matrixStats::rowMaxs(X, na.rm = TRUE)
      } else {
        m <- apply(X, 1L, max, na.rm = TRUE)
      }
      m[!is.finite(m) | m <= 0] <- 1
      if (isTRUE(sparse)) {
        X <- Matrix::Matrix(X, sparse = TRUE)
        X <- Matrix::Diagonal(x = 1 / m) %*% X
        if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
      } else {
        X <- sweep(X, 1L, m, "/")
      }
    } else if (isTRUE(sparse)) {
      X <- Matrix::Matrix(X, sparse = TRUE)
      if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
    }

  } else if (is.data.frame(peaks)) {
    X <- as.matrix(peaks)
    storage.mode(X) <- "double"
    if (normalize == "max") {
      if (requireNamespace("matrixStats", quietly = TRUE)) {
        m <- matrixStats::rowMaxs(X, na.rm = TRUE)
      } else {
        m <- apply(X, 1L, max, na.rm = TRUE)
      }
      m[!is.finite(m) | m <= 0] <- 1
      if (isTRUE(sparse)) {
        X <- Matrix::Matrix(X, sparse = TRUE)
        X <- Matrix::Diagonal(x = 1 / m) %*% X
        if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
      } else {
        X <- sweep(X, 1L, m, "/")
      }
    } else if (isTRUE(sparse)) {
      X <- Matrix::Matrix(X, sparse = TRUE)
      if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
    }

  } else if (is.list(peaks)) {
    # list of numeric vectors (aligned by names if present)
    if (length(peaks) == 0L) stop("peaks list is empty")

    # pre-normalize vectors if requested (saves a big matrix pass)
    if (normalize == "max") {
      m <- vapply(peaks, function(v) {
        v <- as.numeric(v)
        if (length(v)) max(v, na.rm = TRUE) else 0
      }, numeric(1))
      m[!is.finite(m) | m <= 0] <- 1
      peaks <- Map(function(v, s) as.numeric(v) / s, peaks, m)
    }

    nms_list <- lapply(peaks, names)
    have_names <- vapply(nms_list, function(x) !is.null(x) && length(x) > 0L, logical(1))

    if (all(have_names)) {
      all_names <- sort(unique(unlist(nms_list, use.names = FALSE)))
      X <- do.call(rbind, lapply(peaks, function(v) {
        out <- numeric(length(all_names))
        names(out) <- all_names
        idx <- match(names(v), all_names)
        out[idx] <- as.numeric(v)
        out
      }))
      colnames(X) <- all_names
    } else {
      k <- unique(vapply(peaks, length, integer(1)))
      if (length(k) != 1L) stop("Unnamed list elements must all have the same length.")
      X <- do.call(rbind, lapply(peaks, as.numeric))
    }

    storage.mode(X) <- "double"
    if (isTRUE(sparse)) {
      X <- Matrix::Matrix(X, sparse = TRUE)
      if (!methods::is(X, "dgCMatrix")) X <- methods::as(X, "dgCMatrix")
    }

  } else {
    stop("peaks must be a list of MALDIquant::MassPeaks, a numeric matrix/data.frame, or a list of numeric vectors.")
  }

  # Label checks
  if (nrow(X) != length(labels)) {
    stop("length(labels) must equal number of samples (rows of X).")
  }

  list(X = X, Y = labels)
}
