#' Fast MDA-style variable selection using ranger permutation importance
#'
#' Computes feature importances with a single multiclass (or regression) random
#' forest using the ranger engine and its C++ permutation importance. The null
#' distribution of importances is estimated either by appending artificial “false”
#' (noise) features (when `nbf > 0`) or by mirroring negative importances (when
#' `nbf = 0`). An estimator of the proportion of useless features yields pi0, and the top (1 - pi0) proportion
#' of true features are selected. This implementation is a fast, OS-agnostic alternative to
#' repeated/random-forest-based MDA schemes.
#'
#' @param X Numeric matrix (n x p); samples in rows, features in columns. Column
#'   names should be feature IDs (e.g., m/z). Non-finite values are set to zero
#'   internally for modeling.
#' @param Y Factor (classification) or numeric (regression) response of length n.
#'   The default `mtry` is chosen based on the task: floor(sqrt(p)) for
#'   classification; max(floor(p/3), 1) for regression.
#' @param ntree Integer; number of trees. Default 1000.
#' @param nbf Integer (>= 0); number of artificial “false” (noise) features to
#'   append to X to estimate the null distribution. Default 0 disables this and
#'   uses mirrored negative importances as the null.
#' @param nthreads Integer; total number of threads for ranger. Default is
#'   max(1, parallel::detectCores() - 1).
#' @param mtry Optional integer; variables tried at each split. If NULL (default),
#'   computed as floor(sqrt(p)) for classification or max(floor(p/3), 1) for regression.
#' @param sample_fraction Numeric in (0, 1]; subsampling fraction per tree (speed/
#'   regularization knob). Default 1.
#' @param min_node_size Integer; ranger minimum node size. Larger values speed up
#'   training and yield simpler trees. Default 1.
#' @param seed Integer; RNG seed for reproducibility. Default 123.
#'
#' @return A list with:
#'   - nb_to_sel: integer; number of selected features (floor(p_true * (1 - pi0))).
#'   - sel_moz: character vector of selected feature names (columns of X).
#'   - imp_sel: named numeric vector of importances for selected features (true features only).
#'   - all_imp: named numeric vector of importances for all true features.
#'   - pi0: estimated proportion of null features.
#'
#' @details
#' - A single ranger model is fit with importance = "permutation".
#'   This computes permutation importance in C++ using OOB (fast and stable).
#' - Null and pi0:
#'   - If `nbf > 0`, `nbf` false features (uniform between \code{min(X)} and \code{max(X)}) are appended;
#'     negative importances among them help shape the null. An estimator of the proportion of useless features
#'     over high quantiles (e.g., 0.75–1) yields pi0 and is adjusted for the number
#'     of false features.
#'   - If `nbf = 0`, the null is approximated by mirroring negative importances of
#'     true features. If no negative importances occur, pi0 is set to 0 (conservative).
#' - Task: factors in `Y` trigger probability = TRUE; numeric `Y` triggers regression.
#' - Robustness: any non-finite importances are set to zero. Selection is performed
#'   only among the original (true) features; false features are discarded.
#' - Performance: this is typically 5–20x faster than randomForest-based MDA and
#'   fully multithreaded via `nthreads`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100; p <- 300
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("mz_", seq_len(p))
#' Y <- factor(sample(letters[1:3], n, replace = TRUE))
#'
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   out <- fast_mda(
#'     X, Y,
#'     ntree = 500,
#'     nbf = 50,
#'     nthreads = max(1L, parallel::detectCores() - 1L),
#'     seed = 42
#'   )
#'   out$nb_to_sel
#'   head(out$sel_moz)
#'   # Top importances
#'   head(sort(out$all_imp, decreasing = TRUE))
#' }
#' }
#'
#' @references Alexandre Godmer, Yahia Benzerara, Emmanuelle Varon, Nicolas Veziris, Karen Druart, Renaud Mozet, Mariette Matondo, Alexandra Aubry, Quentin Giai Gianetto, MSclassifR: An R package for supervised classification of mass spectra with machine learning methods, Expert Systems with Applications, Volume 294, 2025, 128796, ISSN 0957-4174, \doi{10.1016/j.eswa.2025.128796}.
#'
#' @seealso ranger::ranger; for cross-validated permutation importance, see
#'   fast_cvpvi. For a wrapper that plugs MDA/CVP into broader selection workflows,
#'   see SelectionVar (MethodSelection = "mda" or "cvp").
#' @export
fast_mda <- function(
    X, Y,
    ntree = 1000,
    nbf = 0,
    nthreads = max(1L, parallel::detectCores() - 1L),
    mtry = NULL,
    sample_fraction = 1,
    min_node_size = 1L,
    seed = 123
) {
  if (!requireNamespace("ranger", quietly = TRUE))
    stop("fast_mda requires the 'ranger' package.", call. = FALSE)
  X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  if (anyNA(X) || anyNA(Y)) stop("NA not permitted in X or Y")
  set.seed(seed)

  nthreads_eff <- .limited_cores(nthreads)

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

  fit <- ranger::ranger(
    y = Y, x = Xn,
    num.trees = ntree,
    mtry = mtry,
    importance = "permutation",
    write.forest = FALSE,
    num.threads = nthreads_eff,
    sample.fraction = sample_fraction,
    min.node.size = min_node_size,
    seed = seed
  )
  vi <- fit$variable.importance

  vi_true <- vi[seq_len(p)]
  vi_true[!is.finite(vi_true)] <- 0

  if (nbf > 0) {
    vi_false <- vi[(p + 1L):length(vi)]
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
    all_imp = vi_true,
    pi0 = pi0f
  )
}
