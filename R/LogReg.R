#' Fast supervised classifier with m/z subsetting and optional sampling
#'
#' Trains a multiclass classifier on a subset of m/z features using cross-validation.
#' For kind = "rf", it automatically delegates to a ranger-based algorithm
#' (LogReg_rf_fast) when available for maximum speed and parallelism; otherwise it
#' uses the caret R package  with method = "ranger" as a fast fallback. Other kinds ("linear",
#' "nnet", "svm", "xgb") are trained via caret with compact grids and optional
#' parallelization. Features (columns) are selected by matching their numeric
#' column names to `moz`. Optional class-balancing (among up/down-sampling or SMOTE) can be applied.
#'
#' @param X Numeric matrix or data.frame with samples in rows and features (m/z) in columns.
#'   Column names must be numeric (or coercible), e.g., "1234.567" or "mz_1234.567".
#'   Non-finite values are set to 0.
#' @param moz Numeric vector of m/z values to keep. Only columns of `X` whose
#'   numeric names match values in `moz` are used. An error is thrown if none match.
#' @param Y Factor (or coercible) of class labels; length must equal nrow(X).
#' @param number Integer; number of CV folds (k). Default 2.
#' @param repeats Integer; number of CV repeats. Default 2.
#' @param Metric Character; selection metric. One of "Kappa", "Accuracy", "F1",
#'   "AdjRankIndex", "MatthewsCorrelation". For non-caret metrics, custom summary
#'   functions are used.
#' @param kind Character; model type. One of "linear" (multinom), "nnet" (nnet),
#'   "rf" (random forest), "svm" (svmLinear2), "xgb" (xgbTree). Default "linear".
#' @param Sampling Character; class-balancing strategy. One of "no", "up", "down",
#'   "smote". For "smote", the function smote_classif(Y ~ ., data.frame(Y, X))
#'   is used before training. For "up"/"down", caret’s in-fold sampling is used.
#' @param ncores Integer; number of CPU cores to use for caret’s parallel backend
#'   (doParallel). Default is all but one core. Ignored if doParallel is unavailable.
#' @param num.trees Integer; number of trees for random forests (ranger engine). Default 500.
#'   Used when kind = "rf" and either the caret "ranger" fallback is used or
#'   the caret-free LogReg_rf_fast is available.
#' @param tuneLength Integer; size of the hyperparameter search (caret-based models).
#'   Default 5 (compact grid).
#' @param seed Integer; random seed for reproducibility. Default 123.
#'
#' @return A list with:
#'   - train_mod: the fitted model (caret::train object) or, if kind = "rf" and
#'     LogReg_rf_fast is available, the structure returned by LogReg_rf_fast
#'     (contains the final ranger model and CV details).
#'   - boxplot: ggplot object of resampling metric distributions (caret paths) or
#'     the boxplot returned by LogReg_rf_fast.
#'   - Confusion.Matrix: caret::confusionMatrix on the fitted model (caret paths) or
#'     the confusion matrix returned by LogReg_rf_fast.
#'   - stats_global: data.frame summarizing per-fold metrics (Metric, Mean, Sd) for
#'     caret paths; from LogReg_rf_fast otherwise.
#'
#' @details
#' - Feature subsetting: `X` is subset to columns whose numeric names match `moz`.
#'   This avoids expensive joins/transposes and guarantees stable feature order.
#' - Random forests: if the function LogReg_rf_fast is available in the namespace
#'   (see its documentation), this function delegates the "rf" case to it for
#'   maximum speed and Windows-friendly parallel CV. Otherwise, it uses caret with
#'   method = "ranger" (still fast and parallelizable).
#' - Sampling: "smote" is applied once, before training; "up"/"down" are applied
#'   in-fold by caret via trainControl(sampling = ...). "no" leaves the data unchanged.
#' - Parallelism: if ncores > 1 and doParallel is installed, a PSOCK cluster is
#'   registered for caret. The fast RF engine (LogReg_rf_fast) internally handles
#'   fold-level parallelism and ranger threading to avoid oversubscription.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(runif(2000), nrow = 100, ncol = 20)
#' colnames(X) <- as.character(round(seq(1000, 1190, length.out = 20), 4))
#' moz <- as.numeric(colnames(X))[seq(1, 20, by = 2)]
#' Y <- factor(sample(letters[1:3], 100, replace = TRUE))
#'
#' # Fast RF (delegates to LogReg_rf_fast if available; else caret + ranger)
#' fit_rf <- LogReg(X, moz, Y, number = 3, repeats = 1, kind = "rf",
#'                  Metric = "Kappa", Sampling = "no", ncores = 4,
#'                  num.trees = 300, seed = 42)
#' fit_rf$Confusion.Matrix
#'
#' # Linear (multinom) with macro F1 metric
#' fit_lin <- LogReg(X, moz, Y, number = 3, repeats = 1, kind = "linear",
#'                   Metric = "F1", Sampling = "no", ncores = 2)
#' fit_lin$stats_global
#' }
#'
#' @seealso LogReg_rf_fast, ranger::ranger, caret::train, caret::confusionMatrix
#' @export
#'
LogReg <- function(X,
                   moz,
                   Y,
                   number = 2,
                   repeats = 2,
                   Metric = c("Kappa", "Accuracy", "F1", "AdjRankIndex", "MatthewsCorrelation"),
                   kind = "linear",
                   Sampling = c("no", "up", "down", "smote"),
                   ncores = max(1L, parallel::detectCores() - 1L),
                   num.trees = 500L,
                   tuneLength = 5L,
                   seed = 123L) {

  if (Metric == "AdjRankIndex" && !requireNamespace("mclust", quietly = TRUE)){
    stop("Metric 'AdjRankIndex' requires 'mclust'.");}
  if (Metric == "MatthewsCorrelation" && !requireNamespace("mltools", quietly = TRUE)){
    stop("Metric 'MatthewsCorrelation' requires 'mltools'.");}

  message("LogReg function according to the following parameters:")
  set.seed(seed)
  
  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "") {
    ncores <- 1L
  }
  ncores <- as.integer(max(1L, ncores))

  if (inherits(X, "Matrix")) X <- as.matrix(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0

  moz <- as.numeric(moz)
  xcols_num <- suppressWarnings(as.numeric(colnames(X)))
  if (anyNA(xcols_num)) {
    xcols_num <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", colnames(X))))
  }
  if (anyNA(xcols_num)) {
    stop("Column names of X must be numeric m/z (or coercible after stripping non-numeric chars).")
  }

  idx <- match(moz, xcols_num)
  keep <- which(!is.na(idx))
  if (length(keep) == 0L) stop("None of the provided moz were found in colnames(X).")

  Data_cross <- X[, idx[keep], drop = FALSE]
  colnames(Data_cross) <- as.character(moz[keep])
  rownames(Data_cross) <- rownames(X)
  if (anyDuplicated(colnames(Data_cross))) {
    dup <- duplicated(colnames(Data_cross))
    Data_cross <- Data_cross[, !dup, drop = FALSE]
    warning("Duplicate m/z columns found after subsetting; keeping first occurrences.")
  }

  Y_target <- factor(Y, labels = make.names(levels(factor(Y))))
  DFnnet <- data.frame(Y_target = Y_target, Data_cross, check.names = FALSE)

  f1_score <- function(predicted, expected, positive.class = names(which.min(table(expected)))) {
    expected  <- factor(expected)
    predicted <- factor(as.character(predicted), levels = levels(expected))
    cm <- as.matrix(table(expected, predicted))
    if (nrow(cm) == 0 || ncol(cm) == 0) return(0)
    precision <- diag(cm) / pmax(1, colSums(cm))
    recall    <- diag(cm) / pmax(1, rowSums(cm))
    f1c <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    f1c[is.na(f1c)] <- 0
    if (nlevels(expected) == 2) f1c[positive.class] else mean(f1c)
  }
  f1 <- function(data, lev = NULL, model = NULL) { val <- f1_score(data$pred, data$obs); names(val) <- "F1"; val }
  adjRankIndex <- function(data, lev = NULL, model = NULL) {
    ARI_val <- mclust::adjustedRandIndex(data$pred, data$obs); names(ARI_val) <- "AdjRankIndex"; ARI_val
  }
  matthewsCorrelation <- function(data, lev = NULL, model = NULL) {
    Mc_val <- mltools::mcc(data$pred, data$obs); names(Mc_val) <- "MatthewsCorrelation"; Mc_val
  }

  Sampling <- match.arg(Sampling)
  Metric   <- match.arg(Metric)
  kind     <- match.arg(kind, c("linear", "nnet", "rf", "svm", "xgb"))

  sampling_for_caret <- NULL
  if (kind == "rf") {
    message(if (Sampling == "no") "No sampling method selected" else paste("Sampling method:", Sampling))
  } else {
    if (Sampling == "smote") {
      warning("SMOTE before CV can inflate metrics for non-RF caret models. ",
              "Using in-fold up-sampling instead to avoid leakage. ",
              "For leakage-free SMOTE, prefer kind='rf' (fast engine).")
      sampling_for_caret <- "up"
    } else if (Sampling %in% c("up", "down")) {
      sampling_for_caret <- Sampling
    } else {
      sampling_for_caret <- NULL
    }
    message(if (is.null(sampling_for_caret)) "No sampling method selected"
            else paste("Sampling method:", sampling_for_caret, "(in-fold)"))
  }

  use_summary <- if (Metric %in% c("Kappa", "Accuracy")) {
    caret::defaultSummary
  } else if (Metric == "F1") {
    f1
  } else if (Metric == "AdjRankIndex") {
    adjRankIndex
  } else if (Metric == "MatthewsCorrelation") {
    matthewsCorrelation
  }

  cores <- .limited_cores(ncores)
  allow_parallel <- FALSE
  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "") {
    allow_parallel <- cores > 1L
  }

  fit.control <- caret::trainControl(
    method = "repeatedcv",
    number = number,
    repeats = repeats,
    classProbs = TRUE,
    summaryFunction = use_summary,
    allowParallel = allow_parallel,
    sampling = sampling_for_caret
  )

  if (allow_parallel &&
      Sys.getenv("_R_CHECK_LIMIT_CORES_", "") == "" &&
      requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    if (requireNamespace("foreach", quietly = TRUE)) {
      on.exit({ try(parallel::stopCluster(cl), silent = TRUE); foreach::registerDoSEQ() }, add = TRUE)
    } else {
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    }
  } else if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::registerDoSEQ()
  }

  if (kind == "rf") {
    if (requireNamespace("ranger", quietly = TRUE) && exists("LogReg_rf_fast", mode = "function")) {
      message("Using fast ranger engine (caret-free) for RF")
      res <- LogReg_rf_fast(
        X = X, moz = moz,
        Y = Y,
        number = number, repeats = repeats,
        Metric = Metric, Sampling = Sampling,
        ncores = if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "") 1L else cores,
        num.trees = num.trees,
        tuneLength = tuneLength, seed = seed
      )
      return(res)
    } else {
      message("Using caret + ranger (fallback)")
      p <- ncol(DFnnet) - 1L
      default_mtry <- floor(sqrt(max(1L, p)))
      tuneGrid <- expand.grid(
        mtry = unique(pmax(1L, round(seq(default_mtry * 0.5, default_mtry * 1.5, length.out = max(2L, tuneLength))))),
        splitrule = "gini",
        min.node.size = c(1L, 5L, 10L)[seq_len(min(3L, tuneLength))]
      )
      modelCV <- caret::train(
        Y_target ~ ., data = DFnnet,
        method = "ranger",
        trControl = fit.control,
        metric = Metric,
        preProcess = NULL,
        tuneGrid = tuneGrid,
        num.trees = as.integer(num.trees),
        num.threads = if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "") 1L else as.integer(cores),
        importance = "none",
        probability = TRUE,
        respect.unordered.factors = "order"
      )
    }

  } else if (kind == "linear") {
    message("Estimation with linear (multinom) method")
    modelCV <- caret::train(
      Y_target ~ ., data = DFnnet, method = "multinom",
      trControl = fit.control, metric = Metric,
      preProcess = c("center", "scale"),
      trace = FALSE, maxit = 1000, MaxNWts = 84581
    )

  } else if (kind == "nnet") {
    message("Estimation with nnet method")
    modelCV <- caret::train(
      Y_target ~ ., data = DFnnet, method = "nnet",
      trControl = fit.control, metric = Metric,
      preProcess = c("center", "scale"),
      trace = FALSE, maxit = 1000, MaxNWts = 84581,
      tuneLength = max(3L, tuneLength)
    )

  } else if (kind == "svm") {
    message("Estimation with SVM (linear)")
    modelCV <- caret::train(
      Y_target ~ ., data = DFnnet, method = "svmLinear2",
      trControl = fit.control, metric = Metric,
      preProcess = c("center", "scale"),
      tuneLength = max(3L, tuneLength)
    )

  } else if (kind == "xgb") {
    message("Estimation with xgbTree method (compact grid)")
    
    # FIX R-DEVEL : forcer la matérialisation (évite erreur ALTREP)
    DFnnet <- as.data.frame(DFnnet)
    names(DFnnet) <- paste0(names(DFnnet), "")
    
    xgbGrid <- expand.grid(
      nrounds = c(25, 50, 100),
      max_depth = c(3, 6),
      eta = c(0.1, 0.3),
      gamma = 0,
      colsample_bytree = c(0.7, 1),
      min_child_weight = c(1, 5),
      subsample = c(0.8, 1)
    )
    
    modelCV <- caret::train(
      Y_target ~ ., data = DFnnet, method = "xgbTree",
      trControl = fit.control, metric = Metric,
      preProcess = NULL, tuneGrid = xgbGrid, verbose = FALSE
    )
    
  } else {
    stop("Unknown kind: ", kind)
  }

  variable <- value <- NULL
  b_1 <- suppressMessages(reshape2::melt(modelCV$resample))
  b1 <- ggplot2::ggplot(data = b_1, ggplot2::aes(variable, value, color = variable)) +
    ggplot2::geom_boxplot() + ggplot2::theme_bw() + ggplot2::ylim(0, 1)

  Mean_Metric <- dplyr::summarize(dplyr::group_by(b_1, variable), Mean = mean(value, na.rm = TRUE))
  Sd_metric   <- dplyr::summarize(dplyr::group_by(b_1, variable), Sd = stats::sd(value, na.rm = TRUE))
  statsGlobal <- merge(Mean_Metric, Sd_metric)
  colnames(statsGlobal) <- c("Metric", "Mean", "Sd")

  Confusion.Matrix <- caret::confusionMatrix(modelCV, "none")

  list(
    train_mod = modelCV,
    boxplot = b1,
    Confusion.Matrix = Confusion.Matrix,
    stats_global = statsGlobal
  )
}
