#' Fast random-forest classifier with stratified CV and in-fold sampling (ranger, caret-free)
#'
#' Trains a multiclass random-forest classifier using the ranger algorithm with a
#' compact hyperparameter search and repeated stratified cross-validation. Feature
#' columns are first subset by the provided m/z list (moz). Class balancing
#' (no/up/down/SMOTE) is applied only within training folds to avoid leakage, and
#' again on the full data before fitting the final model. The evaluation across
#' folds can be parallelized in a Windows-safe manner (PSOCK), while avoiding
#' CPU oversubscription by giving each fold worker one ranger thread. Returns the
#' final ranger model, per-fold metrics, a confusion matrix on the full data, and
#' a ggplot boxplot of resampling metrics.
#'
#' @param X Numeric matrix or data frame; rows are samples and columns are features (m/z).
#'   Column names must be numeric (coercible with as.numeric), representing the feature m/z.
#'   Non-finite values are set to 0 internally.
#' @param moz Numeric vector of m/z to keep. Only columns of X whose numeric names
#'   match values in `moz` are used. An error is raised if none match.
#' @param Y Factor (or coercible to factor) of class labels; length must equal nrow(X).
#' @param number Integer; number of CV folds (k). Default 5.
#' @param repeats Integer; number of CV repeats. Default 1.
#' @param Metric Character; CV selection metric. One of
#'   "Kappa", "Accuracy", "F1", "AdjRankIndex", "MatthewsCorrelation".
#'   The best hyperparameters maximize this metric averaged over folds.
#' @param Sampling Character; class-balancing strategy applied within each training fold
#'   (and before the final fit on the full data). One of "no", "up", "down", "smote".
#'   - "up": up-samples minority classes to the majority count (base R).
#'   - "down": down-samples majority classes to the minority count (base R).
#'   - "smote": uses the package’s internal smote_classif(Y ~ ., data.frame(Y, X), C.perc = "balance").
#' @param ncores Integer; number of CPU cores to use. Controls both fold-level parallelism
#'   and ranger threads when not parallelizing folds. Default is all but one core.
#' @param num.trees Integer; number of trees per ranger model. Default 500.
#' @param tuneLength Integer; upper bound on the size of the hyperparameter grid.
#'   If the full grid (mtry × min.node.size) is larger, a random subset of size
#'   `tuneLength` is used. Default 5.
#' @param folds_parallel Character; "auto", "TRUE", or "FALSE".
#'   - "auto": parallelize across folds when ncores >= 2 and total folds (number × repeats) >= 2.
#'   - "TRUE": force fold-level parallelism (PSOCK on Windows).
#'   - "FALSE": evaluate folds sequentially; ranger then uses up to `ncores` threads per fit.
#' @param seed Integer; RNG seed for reproducibility. Default 123.
#' @param mtry Optional integer; if provided, fixes the number of variables tried at each split.
#'   If NULL (default), a small grid around floor(sqrt(p)) is used, where p = number of features.
#' @param splitrule Character; ranger split rule (e.g., "gini", "extratrees"). Default "gini".
#' @param sample.fraction Numeric in (0, 1]; subsampling fraction per tree in ranger. Default 1.
#' @param min.node.size.grid Integer vector; candidate values for ranger’s `min.node.size`
#'   used to build the tuning grid. Default c(1, 5, 10).
#' @param min_node_frac Numeric in (0, 1]. Safety cap for ranger’s min.node.size
#'   per fold/final fit: the value used is min(requested_min.node.size,
#'   floor(min_node_frac * n_train)), with a lower bound of 1.
#'   This prevents root-only trees (near-uniform class probabilities) on small
#'   training folds (e.g., with SMOTE). Applied inside CV and for the final model.
#'   Default: 1/3 (set to 1 to disable capping).
#'
#' @return A list with:
#'   - train_mod: list with fields
#'       - model: the fitted ranger::ranger object (final model on full data)
#'       - method: "ranger"
#'       - best_params: data.frame with the best hyperparameters found by CV
#'       - cv_score: best mean CV score (according to `Metric`)
#'       - metric: the metric name used
#'   - boxplot: ggplot object showing the distribution of per-fold metric values
#'   - Confusion.Matrix: caret::confusionMatrix for predictions of the final model on the full data
#'   - stats_global: data.frame with columns Metric, Mean, Sd summarizing per-fold metrics
#'   - resample: data.frame of per-fold metrics (columns: variable, value, fold)
#'
#' @details
#' - Feature subsetting: X is subset to columns whose numeric names match `moz`. This avoids
#'   expensive joins/transposes and guarantees consistent feature order.
#' - Cross-validation: folds are stratified by Y and repeated `repeats` times. Sampling is applied
#'   only to training indices in each fold (to prevent leakage) and again before the final fit.
#' - Hyperparameter search: a compact grid over mtry (around sqrt(p)) and min.node.size
#'   (from `min.node.size.grid`), optionally downsampled to `tuneLength`. The best combination
#'   maximizes the chosen metric averaged over folds.
#' - Parallel strategy: by default ("auto"), the code parallelizes across folds with a PSOCK cluster
#'   (Windows-safe) and sets ranger’s num.threads = 1 inside each worker to avoid oversubscription.
#'   If you set folds_parallel = "FALSE", folds run sequentially and each ranger fit uses up to
#'   `ncores` threads for strong single-fit parallelism.
#' - Metrics:
#'   - Accuracy and Cohen’s Kappa computed from the confusion matrix.
#'   - F1 is macro-averaged across classes.
#'   - AdjRankIndex uses mclust::adjustedRandIndex.
#'   - MatthewsCorrelation is the multiclass MCC.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(runif(3000), nrow = 100, ncol = 30)
#' colnames(X) <- as.character(round(seq(1000, 1290, length.out = 30), 4))
#' moz <- as.numeric(colnames(X))[seq(1, 30, by = 2)]  # keep half the m/z
#' Y <- factor(sample(letters[1:3], 100, replace = TRUE))
#'
#' fit <- LogReg_rf_fast(
#'   X, moz, Y,
#'   number = 3, repeats = 1,
#'   Metric = "Kappa",
#'   Sampling = "no",
#'   ncores = 4,
#'   num.trees = 300,
#'   tuneLength = 4,
#'   seed = 42
#' )
#' fit$train_mod$best_params
#' fit$Confusion.Matrix
#' }
#'
#' @seealso ranger::ranger, caret::confusionMatrix
#' @export
LogReg_rf_fast <- function(
    X,
    moz,
    Y,
    number = 5,
    repeats = 1,
    Metric = c("Kappa", "Accuracy", "F1", "AdjRankIndex", "MatthewsCorrelation"),
    Sampling = c("no", "up", "down", "smote"),
    ncores = max(1L, parallel::detectCores() - 1L),
    num.trees = 500L,
    tuneLength = 5L,
    folds_parallel = c("auto", "TRUE", "FALSE"),
    seed = 123L,
    mtry = NULL,
    splitrule = "gini",
    sample.fraction = 1.0,
    min.node.size.grid = c(1L, 5L, 10L),
    min_node_frac = 1/3  # cap min.node.size to at most floor(min_node_frac * n_train)
) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("LogReg_rf_fast requires the 'ranger' package. Please install it.", call. = FALSE)
  }
  set.seed(seed)
  Metric <- match.arg(Metric)
  Sampling <- match.arg(Sampling)
  folds_parallel <- match.arg(as.character(folds_parallel), c("auto", "TRUE", "FALSE"))

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  moz <- as.numeric(moz)
  xcols_num <- suppressWarnings(as.numeric(colnames(X)))
  if (anyNA(xcols_num)) stop("Column names of X must be numeric (m/z).")
  idx <- match(moz, xcols_num)
  keep <- which(!is.na(idx))
  if (length(keep) == 0L) stop("None of the provided moz were found in colnames(X).")
  X <- X[, idx[keep], drop = FALSE]
  colnames(X) <- as.character(moz[keep])
  X[!is.finite(X)] <- 0
  Y <- factor(Y)
  if (nrow(X) != length(Y)) stop("nrow(X) must match length(Y).")

  n <- nrow(X); p <- ncol(X)
  levelsY <- levels(Y)

  create_folds <- function(Y, k, repeats, seed) {
    set.seed(seed)
    nn <- seq_along(Y)
    folds_out <- vector("list", 0)
    if (requireNamespace("caret", quietly = TRUE)) {
      for (r in seq_len(repeats)) {
        fl <- caret::createFolds(Y, k = k, list = TRUE, returnTrain = FALSE)
        for (f in seq_len(k)) {
          test_idx <- fl[[f]]; if (length(test_idx) == 0L) next
          train_idx <- setdiff(nn, test_idx)
          folds_out[[length(folds_out) + 1L]] <- list(train = train_idx, test = test_idx)
        }
      }
    } else {
      by_class <- split(nn, Y)
      for (r in seq_len(repeats)) {
        ixc <- lapply(by_class, function(ix) sample(ix))
        folds_class <- lapply(ixc, function(ix) {
          if (length(ix) == 0L) vector("list", k) else {
            sp <- split(ix, rep_len(seq_len(k), length(ix)))
            if (length(sp) < k) sp[(length(sp) + 1L):k] <- list(integer(0)); sp
          }
        })
        for (f in seq_len(k)) {
          test_idx <- unlist(lapply(folds_class, function(fc) fc[[f]]), use.names = FALSE)
          if (length(test_idx) == 0L) next
          train_idx <- setdiff(nn, test_idx)
          folds_out[[length(folds_out) + 1L]] <- list(train = train_idx, test = test_idx)
        }
      }
    }
    if (length(folds_out) == 0L) stop("Failed to create CV folds; reduce 'number' or 'repeats'.")
    folds_out
  }
  folds <- create_folds(Y, number, repeats, seed)
  K <- length(folds)

  upsample_matrix <- function(X, Y) {
    tab <- table(Y); maxn <- max(tab)
    idx_list <- split(seq_along(Y), Y)
    idx_new <- unlist(lapply(idx_list, function(ix) {
      if (length(ix) < maxn) c(ix, sample(ix, maxn - length(ix), replace = TRUE)) else ix
    }), use.names = FALSE)
    idx_new <- sample(idx_new, length(idx_new))
    list(X = X[idx_new, , drop = FALSE], Y = factor(Y[idx_new], levels = levelsY))
  }
  downsample_matrix <- function(X, Y) {
    tab <- table(Y); minn <- min(tab)
    idx_list <- split(seq_along(Y), Y)
    idx_new <- unlist(lapply(idx_list, function(ix) sample(ix, minn)), use.names = FALSE)
    idx_new <- sample(idx_new, length(idx_new))
    list(X = X[idx_new, , drop = FALSE], Y = factor(Y[idx_new], levels = levelsY))
  }
  smote_df <- function(X, Y) {
    mozv <- colnames(X)
    df <- data.frame(Y = Y, X, check.names = FALSE)
    Smoted <- smote_classif(Y ~ ., df, C.perc = "balance")
    XX <- as.matrix(Smoted[, -1, drop = FALSE]); XX[!is.finite(XX)] <- 0; colnames(XX) <- mozv
    YY <- factor(Smoted$Y, levels = levelsY)
    list(X = XX, Y = YY)
  }

  cm_metrics <- function(y_true, y_pred) {
    y_true <- factor(y_true, levels = levelsY)
    y_pred <- factor(y_pred, levels = levelsY)
    cm <- table(y_true, y_pred)
    N <- sum(cm); if (N == 0) N <- 1
    acc <- sum(diag(cm)) / N
    p0 <- acc; rows <- rowSums(cm); cols <- colSums(cm)
    pe <- sum(rows * cols) / (N * N)
    kap <- if (pe < 1) (p0 - pe) / (1 - pe) else 0
    prec <- diag(cm) / pmax(1, cols)
    rec  <- diag(cm) / pmax(1, rows)
    f1c  <- ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
    f1m  <- mean(f1c, na.rm = TRUE)
    ari <- tryCatch(mclust::adjustedRandIndex(y_pred, y_true), error = function(e) NA_real_)
    cdiag <- sum(diag(cm))
    numer <- cdiag * N - sum(rows * cols)
    denom <- sqrt((N^2 - sum(cols^2)) * (N^2 - sum(rows^2)))
    mcc <- if (denom > 0) numer / denom else 0
    c(Accuracy = acc, Kappa = kap, F1 = f1m, AdjRankIndex = ari, MatthewsCorrelation = mcc)
  }

  if (is.null(mtry)) {
    base_mtry <- max(1L, floor(sqrt(p)))
    mtry_grid <- unique(pmax(1L, c(round(base_mtry * 0.5), base_mtry, round(base_mtry * 1.5))))
  } else {
    mtry_grid <- as.integer(pmax(1L, pmin(p, mtry)))
  }
  min_node_grid <- unique(as.integer(pmax(1L, min.node.size.grid)))
  grid <- expand.grid(mtry = mtry_grid, min.node.size = min_node_grid, splitrule = splitrule,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (nrow(grid) > tuneLength) {
    grid <- grid[sample.int(nrow(grid), tuneLength), , drop = FALSE]
  }

  cores <- .limited_cores(ncores)
  do_par <- if (folds_parallel == "auto") (cores >= 2L && K >= 2L) else (folds_parallel == "TRUE")
  fold_workers <- if (do_par) min(cores, K) else 1L
  fold_threads <- if (do_par) 1L else max(1L, cores)

  eval_param <- function(params, keep_predictions = FALSE) {
    set.seed(seed + 17L)
    fold_worker <- function(fidx) {
      tr <- folds[[fidx]]$train; te <- folds[[fidx]]$test
      if (length(te) == 0L) return(list(metric = NA_real_, all_metrics = rep(NA_real_, 5), pred = NULL, te = te))
      Xtr <- X[tr, , drop = FALSE]; Ytr <- Y[tr]
      if (Sampling == "up") {
        tmp <- upsample_matrix(Xtr, Ytr); Xtr <- tmp$X; Ytr <- tmp$Y
      } else if (Sampling == "down") {
        tmp <- downsample_matrix(Xtr, Ytr); Xtr <- tmp$X; Ytr <- tmp$Y
      } else if (Sampling == "smote") {
        tmp <- smote_df(Xtr, Ytr); Xtr <- tmp$X; Ytr <- tmp$Y
      }
      if (length(unique(Ytr)) < 2L) {
        return(list(metric = NA_real_, all_metrics = rep(NA_real_, 5), pred = NULL, te = te))
      }

      min_node_cap <- max(1L, floor(nrow(Xtr) * min_node_frac))
      min_node_use <- max(1L, min(as.integer(params$min.node.size), min_node_cap))

      fit <- ranger::ranger(
        y = Ytr, x = as.data.frame(Xtr, check.names = FALSE),
        num.trees = as.integer(num.trees),
        mtry = as.integer(params$mtry),
        min.node.size = min_node_use,
        splitrule = params$splitrule,
        sample.fraction = sample.fraction,
        probability = TRUE,
        respect.unordered.factors = "order",
        num.threads = as.integer(fold_threads),
        seed = seed + fidx
      )
      pr <- predict(fit, data = as.data.frame(X[te, , drop = FALSE], check.names = FALSE))$predictions
      pred_class <- colnames(pr)[max.col(pr, ties.method = "first")]
      m <- cm_metrics(Y[te], factor(pred_class, levels = levelsY))
      list(metric = m[Metric], all_metrics = m, pred = if (keep_predictions) pred_class else NULL, te = te)
    }

    fold_indices <- seq_len(K)
    fold_results <- if (do_par) {
      cl <- parallel::makeCluster(fold_workers); on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::parLapply(cl, fold_indices, fold_worker)
    } else {
      lapply(fold_indices, fold_worker)
    }

    metric_vals <- vapply(fold_results, function(fr) as.numeric(fr$metric), numeric(1))
    list(score = mean(metric_vals, na.rm = TRUE),
         metric_vals = metric_vals,
         fold_results = fold_results)
  }

  best_score <- -Inf
  best_params <- NULL
  best_eval <- NULL
  for (g in seq_len(nrow(grid))) {
    params <- grid[g, , drop = FALSE]
    ev <- eval_param(params)
    if (!is.finite(ev$score)) next
    if (ev$score > best_score) {
      best_score <- ev$score
      best_params <- params
      best_eval <- ev
    }
  }
  if (is.null(best_params)) stop("Failed to train any model on folds (check sampling or class balance).")

  best_eval <- eval_param(best_params, keep_predictions = TRUE)

  Xfull <- X; Yfull <- Y
  if (Sampling == "up") {
    tmp <- upsample_matrix(Xfull, Yfull); Xfull <- tmp$X; Yfull <- tmp$Y
  } else if (Sampling == "down") {
    tmp <- downsample_matrix(Xfull, Yfull); Xfull <- tmp$X; Yfull <- tmp$Y
  } else if (Sampling == "smote") {
    tmp <- smote_df(Xfull, Yfull); Xfull <- tmp$X; Yfull <- tmp$Y
  }

  min_node_cap_full <- max(1L, floor(nrow(Xfull) * min_node_frac))
  min_node_use_full <- max(1L, min(as.integer(best_params$min.node.size), min_node_cap_full))

  final_model <- ranger::ranger(
    y = Yfull, x = as.data.frame(Xfull, check.names = FALSE),
    num.trees = as.integer(num.trees),
    mtry = as.integer(best_params$mtry),
    min.node.size = min_node_use_full,
    splitrule = best_params$splitrule,
    sample.fraction = sample.fraction,
    probability = TRUE,
    respect.unordered.factors = "order",
    num.threads = as.integer(max(1L, cores)),
    seed = seed
  )

  pr_full <- predict(final_model, data = as.data.frame(X, check.names = FALSE))$predictions
  pred_full <- factor(colnames(pr_full)[max.col(pr_full, ties.method = "first")], levels = levelsY)
  Confusion.Matrix <- caret::confusionMatrix(pred_full, Y)

  resample_df <- do.call(rbind, lapply(seq_len(K), function(fidx) {
    allm <- best_eval$fold_results[[fidx]]$all_metrics
    data.frame(variable = names(allm), value = as.numeric(allm), fold = fidx, stringsAsFactors = FALSE)
  }))

  variable <- value <- NULL
  b1 <- ggplot2::ggplot(resample_df, ggplot2::aes(variable, value, color = variable)) +
    ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::ylim(0, 1)

  statsGlobal <- dplyr::summarize(dplyr::group_by(resample_df, variable),
                                  Mean = mean(value, na.rm = TRUE),
                                  Sd = stats::sd(value, na.rm = TRUE))
  colnames(statsGlobal) <- c("Metric", "Mean", "Sd")

  list(
    train_mod = list(
      model = final_model,
      method = "ranger",
      best_params = best_params,
      cv_score = best_score,
      metric = Metric
    ),
    boxplot = b1,
    Confusion.Matrix = Confusion.Matrix,
    stats_global = statsGlobal,
    resample = resample_df
  )
}
