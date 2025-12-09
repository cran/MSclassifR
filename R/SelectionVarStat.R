#' Fast feature (m/z) selection using multiple hypothesis testing (LIMMA/ANOVA/Kruskal) with optional class balancing (no/up/down/SMOTE)
#'
#' Runs univariate statistical tests across all features (columns) to identify
#' discriminant m/z values. Supports three tests:
#' - "Limma": moderated F-test via the limma R package on k-1 contrasts between groups
#' - "anova": classical one-way ANOVA (implemented in C++ for speed)
#' - "kruskal": Kruskal–Wallis rank-sum test (implemented in C++ with tie correction)
#'
#' Optional class balancing can be applied before testing: no sampling,
#' up-sampling, down-sampling, or SMOTE using the internal \code{smote_classif()} function.
#' P-values are adjusted with the cp4p R package, and features below the FDR
#' threshold are returned, together with an estimated proportion of nulls, \code{pi0}.
#'
#' Missing values are preserved as \code{NA}: for ANOVA/Kruskal, the C++ backends
#' skip NAs per feature; for Limma, missing values are handled by limma internally.
#' Non-finite values (Inf/-Inf) are coerced to \code{NA}.
#'
#' @param X Numeric matrix with samples in rows and features (peaks) in columns.
#'   If a sparse Matrix is provided, it is coerced to a base R dense matrix.
#'   Infinities are set to NA; NAs are allowed.
#' @param Y Class labels. A factor (or coercible to factor) of length \code{nrow(X)}.
#'   Must contain at least two levels.
#' @param stat.test Character string; which univariate test to use. One of
#'   \code{"Limma"}, \code{"anova"}, or \code{"kruskal"}. Default is \code{"Limma"}.
#'   - Limma uses \code{limma::lmFit} on \code{~ 0 + Y}, applies k-1 contrasts
#'     versus a reference group with \code{limma::makeContrasts} and then
#'     \code{limma::eBayes} to compute a moderated F p-value testing overall
#'     between-group differences.
#'   - anova uses a fast C++ one-way ANOVA (upper-tail F p-value).
#'   - kruskal uses a fast C++ Kruskal–Wallis test with tie correction
#'     (upper-tail chi-square p-value).
#' @param pi0.method Character; method for \code{cp4p::estim.pi0} and
#'   \code{cp4p::adjust.p}. Default \code{"abh"}. See cp4p documentation for options.
#' @param fdr Numeric in (0, 1]; FDR threshold used to select features after
#'   p-value adjustment. Default \code{0.05}.
#' @param Sampling Character string; optional class balancing applied before testing.
#'   One of \code{"no"}, \code{"up"}, \code{"down"}, \code{"smote"}. Default \code{"no"}.
#'   - "up": up-samples minority classes to the majority size (base R implementation).
#'   - "down": down-samples majority classes to the minority size (base R).
#'   - "smote": uses the package’s internal \code{smote_classif()} on \code{data.frame(Y, X)}.
#'     Column names of \code{X} are preserved; if SMOTE fails or yields <2 classes,
#'     the function falls back to up-sampling.
#' @param seed Optional integer. If provided, sets the random seed for
#'   up/down sampling and SMOTE to ensure reproducibility.
#'
#' @return A list with:
#'   - \code{nb_to_sel}: integer, \code{floor(ncol(X) * (1 - pi0))} using \code{cp4p::estim.pi0}.
#'   - \code{n_selected_fdr}: integer, number of features with adjusted p-value < \code{fdr}.
#'   - \code{sel_moz}: character vector of selected feature names (columns of \code{X})
#'     with adjusted p-value < \code{fdr}.
#'   - \code{ap}: the full object returned by \code{cp4p::adjust.p} (contains adjusted p-values).
#'
#' @details
#' - Limma design: \code{model.matrix(~ 0 + Y)} is used; a full-rank set of k-1
#'   contrasts vs the first level is constructed with \code{limma::makeContrasts}.
#'   The returned p-values are moderated F-test p-values for overall group differences.
#' - ANOVA/Kruskal implementations are in C++ (see \code{anova_cols_cpp} and
#'   \code{kruskal_cols_cpp}); they process all features in one pass and skip NAs per feature.
#' - Non-finite p-values (if any) are set to 1 before \code{pi0}/adjustment.
#' - SMOTE requires at least two observations per class; otherwise the function
#'   automatically falls back to up-sampling.
#' - For limma, if residual degrees of freedom are not positive
#'   (\code{nrow(X) - nlevels(Y) <= 0}), p-values are set to 1 with a warning.
#'
#' @seealso limma::lmFit, limma::contrasts.fit, limma::eBayes, limma::makeContrasts;
#'   cp4p::estim.pi0, cp4p::adjust.p
#'
#' @examples
#' ###############################################################################
#' ## 1. Pre-processing of mass spectra
#'
#' # load mass spectra and their metadata
#' data("CitrobacterRKIspectra","CitrobacterRKImetadata", package = "MSclassifR")
#' # standard pre-processing of mass spectra
#' spectra <- MSclassifR::SignalProcessing(CitrobacterRKIspectra)
#' # detection of peaks in pre-processed mass spectra
#' peaks <- MSclassifR::PeakDetection(x = spectra, labels = CitrobacterRKImetadata$Strain_name_spot)
#' # build matrix with intensities of peaks (rows = samples, columns = m/z)
#' Y <- factor(CitrobacterRKImetadata$Species)
#' xy <- build_XY_from_peaks(peaks, labels = Y, normalize = "max", sparse = FALSE)
#' X <- xy$X
#' Y <- xy$Y
#'
#' ###############################################################################
#' ## 2. Estimate the optimal number of peaks to discriminate the different species
#'
#' OptiPeaks <- MSclassifR::SelectionVarStat(X,
#'                              Y,
#'                              stat.test = "Limma",
#'                              pi0.method = "abh",
#'                              fdr = 0.05,
#'                              Sampling = "smote",
#'                              seed = 1)
#'
#' ## Estimation of the optimal number of peaks to discriminate species (from the pi0 parameter)
#' OptiPeaks$nb_to_sel
#'
#' ## discriminant mass-to-charge values estimated using a 5 per cent false discovery rate
#' OptiPeaks$sel_moz
#'
#' ## p-values and adjusted p-values estimated for all the tested mass-to-charge values
#' OptiPeaks$ap$adjp
#'
#' @export
SelectionVarStat <- function(X,
                             Y,
                             stat.test = c("Limma", "anova", "kruskal"),
                             pi0.method = "abh",
                             fdr = 0.05,
                             Sampling = c("no", "up", "down", "smote"),
                             seed = NULL) {

  stat.test <- match.arg(stat.test)
  SamplingM <- match.arg(Sampling)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # Coerce inputs
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  Y <- factor(Y)
  stopifnot(nrow(X) == length(Y))
  if (is.null(colnames(X))) colnames(X) <- paste0("mz_", seq_len(ncol(X)))

  # Preserve NAs, only sanitize infinities
  X[is.infinite(X)] <- NA_real_

  if (nlevels(Y) < 2L) stop("Y must contain at least two classes.")

  # Up/down sampling helpers (rows = samples)
  upsample_matrix <- function(X, Y) {
    tab <- table(Y); maxn <- max(tab)
    idx_list <- split(seq_along(Y), Y)
    idx_new <- unlist(lapply(idx_list, function(idx) {
      if (length(idx) < maxn) c(idx, sample(idx, maxn - length(idx), replace = TRUE)) else idx
    }), use.names = FALSE)
    idx_new <- sample(idx_new, length(idx_new))
    list(X = X[idx_new, , drop = FALSE], Y = factor(Y[idx_new], levels = levels(Y)))
  }
  downsample_matrix <- function(X, Y) {
    tab <- table(Y); minn <- min(tab)
    idx_list <- split(seq_along(Y), Y)
    idx_new <- unlist(lapply(idx_list, function(idx) sample(idx, minn, replace = FALSE)),
                      use.names = FALSE)
    idx_new <- sample(idx_new, length(idx_new))
    list(X = X[idx_new, , drop = FALSE], Y = factor(Y[idx_new], levels = levels(Y)))
  }

  # Guard rails for SMOTE
  if (SamplingM == "smote") {
    tab <- table(Y)
    if (length(tab) < 2L || any(tab < 2L)) {
      message("SMOTE requires at least two observations per class; falling back to up-sampling.")
      SamplingM <- "up"
    }
  }

  # Apply optional sampling
  switch(
    SamplingM,
    "no" = {
      message("No sampling method selected")
    },
    "up" = {
      message("Up sampling method selected")
      tmp <- upsample_matrix(X, Y); X <- tmp$X; Y <- tmp$Y
    },
    "down" = {
      message("Down sampling method selected")
      tmp <- downsample_matrix(X, Y); X <- tmp$X; Y <- tmp$Y
    },
    "smote" = {
      message("Smote sampling method selected (internal smote_classif)")
      mozv <- colnames(X)
      df <- data.frame(Y = Y, X, check.names = FALSE)
      Smoted <- try(smote_classif(Y ~ ., df, C.perc = "balance"), silent = TRUE)

      smote_ok <- !(inherits(Smoted, "try-error")) &&
        ("Y" %in% colnames(Smoted)) &&
        length(unique(Smoted$Y)) >= 2L

      if (!smote_ok) {
        message("SMOTE failed or yielded <2 classes; falling back to up-sampling.")
        tmp <- upsample_matrix(X, Y); X <- tmp$X; Y <- tmp$Y
      } else {
        Y <- droplevels(factor(Smoted$Y))
        keep_cols <- mozv[mozv %in% colnames(Smoted)]
        if (length(keep_cols) < length(mozv)) {
          warning("Some features were not found after SMOTE and will be dropped: ",
                  paste(setdiff(mozv, keep_cols), collapse = ", "))
        }
        X <- as.matrix(Smoted[, keep_cols, drop = FALSE])
        storage.mode(X) <- "double"
        colnames(X) <- keep_cols
        X[is.infinite(X)] <- NA_real_
      }
    }
  )

  # Safety checks after sampling
  if (nrow(X) != length(Y)) stop("After sampling, nrow(X) must equal length(Y).")
  Y <- droplevels(factor(Y))
  if (nlevels(Y) < 2L) {
    warning("After sampling there is only one class; returning no selected features.")
    vp <- rep(1, ncol(X))
    ap <- cp4p::adjust.p(vp, pi0.method = pi0.method)
    pi0 <- cp4p::estim.pi0(vp, pi0.method = pi0.method)
    pi0_val <- if (!is.null(pi0$pi0)) as.numeric(pi0$pi0) else as.numeric(pi0[[1]])
    nb_to_sel <- floor(ncol(X) * (1 - pi0_val))
    sel_moz <- character(0L)
    return(list(nb_to_sel = nb_to_sel,
                n_selected_fdr = length(sel_moz),
                sel_moz = sel_moz,
                ap = ap))
  }

  # Compute p-values per feature
  vp <- switch(
    stat.test,
    "Limma" = {
      k <- nlevels(Y)
      if (k < 2L) stop("Need at least 2 classes for Limma.")
      res_df <- nrow(X) - k
      if (res_df <= 0L) {
        warning("Not enough residual degrees of freedom for limma (nrow(X) - nlevels(Y) <= 0). Returning p=1.")
        rep(1, ncol(X))
      } else {
        # Use intercept design; eBayes F-tests the non-intercept coefficients,
        # which is the overall equality of group means.
        design <- stats::model.matrix(~ Y)
        fit <- limma::lmFit(t(X), design)
        fit <- limma::eBayes(fit, robust = TRUE)
        fit$F.p.value
      }
    },
    "anova" = {
      grp <- as.integer(Y) - 1L
      anova_cols_cpp(X, grp)
    },
    "kruskal" = {
      grp <- as.integer(Y) - 1L
      kruskal_cols_cpp(X, grp)
    }
  )

  # Clean p-values
  vp[!is.finite(vp)] <- 1
  vp[vp < 0] <- 0
  vp[vp > 1] <- 1

  # pi0 estimation and FDR adjustment
  pi0 <- cp4p::estim.pi0(vp, pi0.method = pi0.method)
  pi0_val <- if (!is.null(pi0$pi0)) as.numeric(pi0$pi0) else as.numeric(pi0[[1]])
  nb_to_sel <- floor(ncol(X) * (1 - pi0_val))
  ap <- cp4p::adjust.p(vp, pi0.method = pi0.method)

  sel_moz <- colnames(X)[which(ap$adjp$adjusted.p < fdr)]
  list(
    nb_to_sel = nb_to_sel,
    n_selected_fdr = length(sel_moz),
    sel_moz = sel_moz,
    ap = ap
  )
}
