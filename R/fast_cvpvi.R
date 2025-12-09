#' Fast cross-validated permutation variable importance (ranger-based)
#'
#' Computes cross-validated permutation variable importance (PVI) using the
#' ranger random-forest algorithm. For each CV fold, a ranger model is trained on
#' the training split and permutation importance is computed (OOB) inside ranger
#' in C++. Importances are averaged across folds to obtain a stable CV importance
#' vector. Optionally appends artificial “false” features to estimate the null
#' distribution and pi0, then selects the top (1 - pi0) proportion of features.
#' The evaluation can parallelize across folds (Windows-safe PSOCK) while
#' avoiding CPU oversubscription.
#'
#' @param X Numeric matrix (n x p); samples in rows, features in columns. Column
#'   names should be feature IDs (e.g., m/z). Non-finite values are set to zero
#'   internally for modeling.
#' @param Y Factor or numeric response of length n. A factor triggers
#'   classification; numeric triggers regression.
#' @param k Integer; number of cross-validation folds. Default 5.
#' @param ntree Integer; number of trees per fold model. Default 500.
#' @param nbf Integer (>= 0); number of artificial “false” (noise) features to
#'   append to X for estimating the null distribution of importances. Default 0
#'   disables this (the null is then approximated using mirrored negative
#'   importances).
#' @param nthreads Integer; total threads available. When parallelizing folds,
#'   each fold worker gets one ranger thread to avoid oversubscription; when not
#'   parallelizing folds, ranger uses up to `nthreads` threads. Default is
#'   max(1, detectCores() - 1).
#' @param folds_parallel Character; "auto", "TRUE", or "FALSE".
#'   - "auto": parallelize across folds when k > 1 and nthreads >= 4 (default).
#'   - "TRUE": force fold-level parallelism (PSOCK cluster).
#'   - "FALSE": evaluate folds sequentially (ranger can then use multiple threads).
#' @param mtry Optional integer; variables tried at each split. If NULL, defaults
#'   to floor(sqrt(p)) for classification or max(floor(p/3), 1) for regression.
#' @param sample_fraction Numeric in (0, 1]; subsampling fraction per tree (speed/
#'   regularization knob). Default 1.
#' @param min_node_size Integer; ranger minimum node size. Larger values speed up
#'   training and yield smaller trees. Default 1.
#' @param seed Integer; RNG seed. Default 123.
#'
#' @return A list with:
#'   - nb_to_sel: integer; number of selected features (floor(p * (1 - pi0))).
#'   - sel_moz: character vector of selected feature names (columns of X).
#'   - imp_sel: named numeric vector of CV importances for selected features.
#'   - fold_varim: matrix (features x folds) of per-fold permutation importances.
#'   - cv_varim: matrix (features x 1) of averaged importances across folds.
#'   - pi0: estimated proportion of null features.
#'
#' @details
#' - One ranger model is trained per fold (training split). Permutation
#'   importance (importance = "permutation") is computed in C++ using OOB. The
#'   per-fold importances are averaged to obtain CV importances.
#' - Null and pi0: if `nbf > 0`, false peaks are created to get negative importances. For this, `nbf` noise features (uniform between \code{min(X)} and \code{max(X)})
#'   are appended and negative importances among them help shape the null. If
#'   `nbf = 0`, the null is approximated by mirroring negative importances of
#'   true features. An estimator of the proportion of useless features over high quantiles yields pi0.
#'   If no negative importances occur, pi0 is set to 0 (conservative).
#' - Parallelism: with `folds_parallel = "auto"/"TRUE"`, folds run in parallel
#'   using a PSOCK cluster (Windows-safe). Each worker sets ranger num.threads = 1
#'   to avoid oversubscription. With `"FALSE"`, folds are sequential and ranger
#'   uses up to `nthreads` threads, which can be faster for small k or very large p.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 120; p <- 200
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("mz_", seq_len(p))
#' Y <- factor(sample(letters[1:3], n, replace = TRUE))
#'
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   out <- fast_cvpvi(
#'     X, Y,
#'     k = 5,
#'     ntree = 300,
#'     nbf = 50,
#'     nthreads = max(1L, parallel::detectCores() - 1L),
#'     folds_parallel = "auto",
#'     seed = 42
#'   )
#'   head(out$sel_moz)
#'   # CV importances for top features
#'   head(sort(out$cv_varim[,1], decreasing = TRUE))
#' }
#' }
#'
#' @references Alexandre Godmer, Yahia Benzerara, Emmanuelle Varon, Nicolas Veziris, Karen Druart, Renaud Mozet, Mariette Matondo, Alexandra Aubry, Quentin Giai Gianetto, MSclassifR: An R package for supervised classification of mass spectra with machine learning methods, Expert Systems with Applications, Volume 294, 2025, 128796, ISSN 0957-4174, \doi{10.1016/j.eswa.2025.128796}.
#'
#' @seealso ranger::ranger; for a holdout-based (validation-fold) permutation
#'   alternative, see a custom implementation using predict() on permuted
#'   features. For a full feature-selection wrapper, see SelectionVar with
#'   MethodSelection = "cvp".
#'
#' @export
fast_cvpvi <- function(
    X, Y,
    k = 5,
    ntree = 500,
    nbf = 0,
    nthreads = max(1L, parallel::detectCores() - 1L),
    folds_parallel = c("auto", "TRUE", "FALSE"),
    mtry = NULL,
    sample_fraction = 1,
    min_node_size = 1L,
    seed = 123
) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("fast_cvpvi requires the 'ranger' package.", call. = FALSE)
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  if (anyNA(X) || anyNA(Y)) stop("NA not permitted in X or Y")
  set.seed(seed)

  n <- nrow(X); p <- ncol(X)
  is_classif <- is.factor(Y)
  if (is.null(mtry)) mtry <- if (is_classif) floor(sqrt(p)) else max(floor(p/3), 1L)

  if (nbf > 0) {
    X0 <- matrix(stats::runif(nbf * n, min = min(X), max = max(X)), n, nbf)
    colnames(X0) <- paste0("false_", seq_len(nbf))
    Xn <- cbind(X, X0)
  } else {
    Xn <- X
  }

  if (requireNamespace("caret", quietly = TRUE)) {
    folds <- caret::createFolds(Y, k = k, list = TRUE, returnTrain = FALSE)
  } else {
    idx <- sample.int(n); cuts <- round(n / k)
    from <- (0:(k - 1L)) * cuts + 1L
    to <- pmin(seq_len(k) * cuts, n)
    folds <- lapply(seq_len(k), function(i) idx[from[i]:to[i]])
  }

  folds_parallel <- match.arg(as.character(folds_parallel), c("auto", "TRUE", "FALSE"))
  nthreads_eff <- .limited_cores(nthreads)
  do_par <- if (folds_parallel == "auto") (k > 1L && nthreads_eff >= 4L) else (folds_parallel == "TRUE")
  fold_workers <- if (do_par) min(k, nthreads_eff) else 1L
  fold_threads <- if (do_par) 1L else nthreads_eff

  fit_one_fold <- function(i) {
    test_idx <- folds[[i]]
    dY <- Y[-test_idx]
    dX <- Xn[-test_idx, , drop = FALSE]
    fit <- ranger::ranger(
      y = dY, x = dX,
      num.trees = ntree,
      mtry = mtry,
      importance = "permutation",
      write.forest = FALSE,
      num.threads = fold_threads,
      sample.fraction = sample_fraction,
      min.node.size = min_node_size,
      seed = seed + i
    )
    fit$variable.importance
  }

  if (do_par) {
    cl <- parallel::makeCluster(fold_workers); on.exit(parallel::stopCluster(cl), add = TRUE)
    vi_list <- parallel::parLapply(cl, seq_len(k), fit_one_fold)
  } else {
    vi_list <- lapply(seq_len(k), fit_one_fold)
  }

  vi_mat <- do.call(cbind, vi_list)
  colnames(vi_mat) <- paste0(seq_len(k), "-fold")
  rownames(vi_mat) <- colnames(Xn)
  cv_vi <- rowMeans(vi_mat, na.rm = TRUE)

  vi_true <- cv_vi[seq_len(p)]
  vi_true[!is.finite(vi_true)] <- 0

  if (nbf > 0) {
    vi_false <- cv_vi[(p + 1L):length(cv_vi)]
    vi_false_neg <- vi_false[is.finite(vi_false) & vi_false < 0]
    vi1 <- c(vi_true, vi_false_neg)
  } else {
    vi1 <- vi_true
  }

  imp_neg <- vi1[vi1 < 0]
  if (length(imp_neg) == 0L) {
    pi0f <- 0
  } else {
    imp_null <- c(imp_neg, -imp_neg)
    q_ext <- seq(0.75, 1, by = 0.01)
    Fall <- stats::ecdf(vi1)
    pi0_raw <- vapply(q_ext, function(q) {
      qin <- stats::quantile(imp_null, q, na.rm = TRUE)
      min(Fall(qin) / q, 1)
    }, numeric(1))
    if (nbf > 0) {
      Nfn <- sum(vi_false_neg < 0)
      pi0f <- (min(pi0_raw) * (p + Nfn) - Nfn) / p
    } else {
      pi0f <- min(pi0_raw)
    }
  }

  nb_to_sel <- min(p, max(1L, floor(p * (1 - pi0f))))
  sel_idx <- order(-vi_true)[seq_len(nb_to_sel)]
  sel_moz <- names(vi_true)[sel_idx]
  imp_sel <- vi_true[sel_idx]

  list(
    nb_to_sel = nb_to_sel,
    sel_moz = sel_moz,
    imp_sel = imp_sel,
    fold_varim = vi_mat[seq_len(p), , drop = FALSE],
    cv_varim = matrix(vi_true, ncol = 1, dimnames = list(names(vi_true), "CV_PerVarImp")),
    pi0 = pi0f
  )
}
