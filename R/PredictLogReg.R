#' Prediction of the category to which a mass spectrum belongs
#'
#' Predicts the category (species, phenotype, etc.) of each spectrum in a list
#' of MALDIquant MassPeaks using one or more trained models (e.g., multinomial
#' logistic regression from LogReg). Peaks are matched to a given shortlist of
#' discriminant m/z values (moz) within a tolerance; unmatched positions are
#' filled with `noMatch`. If several models are supplied, the function also
#' produces meta-predictions: per-class Fisher combinations (if 'metap' is
#' available) and a majority-vote fraction across models.
#'
#' @param peaks a list of MALDIquant::MassPeaks objects (one per spectrum).
#' @param model a model or a list of models estimated from a shortlist of m/z
#'   (e.g., the output of LogReg). Each model must support predict(..., type = "prob").
#'   If a single model is supplied, it is wrapped into a one-element list.
#' @param moz a numeric vector of shortlisted m/z values used for prediction
#'   (typically the selection from SelectionVar/SelectionVarStat_fast).
#' @param tolerance numeric; accepted m/z tolerance (in Da) for matching peaks to
#'   moz. Default 6.
#' @param toleranceStep numeric; if a spectrum yields no matches at the initial
#'   tolerance, the function retries by increasing tolerance in steps of
#'   `toleranceStep` until at least one match is found (bounded internally). Default 2.
#' @param normalizeFun logical; if TRUE (default), per-spectrum max normalization
#'   is applied after matching (row is divided by its maximum).
#' @param noMatch numeric; intensity used when no peak matches a given m/z. Default 0.
#' @param Reference optional factor of true categories, length equal to length(peaks).
#'   If provided and has at least two distinct levels, the function returns, in
#'   addition to the predictions, per-model confusion matrices (caret::confusionMatrix).
#' @param chunk_size integer; number of spectra per prediction batch (rows of X).
#'   Large datasets can be processed in chunks to limit memory. Default 10000.
#' @param ncores integer; number of cores used while building X from peaks on the
#'   R side (the C++ matching itself is single-threaded here). Default 1.
#' @param verbose logical; print progress messages. Default FALSE.
#'
#' @return If Reference is missing (or has < 2 levels), a data.frame with:
#'   - name: spectrum name (from MassPeaks metaData fullName/file when available)
#'   - method: model identifier (from model$method; suffixed with "_i" if needed)
#'   - one column per class with predicted probabilities
#'   - pred_max_p: predicted class (argmax of probabilities)
#'
#'   If Reference is provided with at least two levels, a list with:
#'   - Prob.results: the predictions data.frame as above
#'   - Confusion.Matrix: a list of caret::confusionMatrix objects (one per method)
#'
#' @details
#' - Matching and normalization: peak-to-moz matching is performed by build_X_from_peaks_fast
#'   (C++-backed; nearest-within-tolerance). If no m/z from a spectrum match the
#'   shortlist initially, the tolerance is increased by `toleranceStep` in a small
#'   number of attempts until at least one match is found. If `normalizeFun = TRUE`,
#'   each row is divided by its maximum (guarded to avoid divide-by-zero).
#' - Multiple models: when several models are supplied, the output contains one
#'   set of probabilities per model (with method column identifying it). Two
#'   additional rows per spectrum can be appended:
#'   - comb_fisher: per-class Fisher combined p-values computed via metap::sumlog (if available).
#'   - max_vote: per-class fraction of models casting the top-probability vote for that class.
#' - Models: this function is agnostic of the modeling engine as long as predict(type = "prob")
#'   is implemented (e.g., caret multinom/nnet/ranger/xgb, glmnet, etc.).
#'
#' @examples
#' \donttest{
#' library(MSclassifR)
#' library(MALDIquant)
#'
#' ## 1) Preprocess and detect peaks
#' data("CitrobacterRKIspectra", "CitrobacterRKImetadata", package = "MSclassifR")
#' spectra <- SignalProcessing(CitrobacterRKIspectra)
#' peaks   <- MSclassifR::PeakDetection(x = spectra, averageMassSpec = FALSE)
#'
#' ## 2) Build X and Y (sample-by-peak intensities + labels)
#' ##    Option A: if you prefer the helper and a sparse return:
#' Y <- factor(CitrobacterRKImetadata$Species)
#' xy <- build_XY_from_peaks(peaks, labels = Y, normalize = "max", sparse = FALSE)
#' X <- xy$X
#' Y <- xy$Y
#'
#' ##    Option B: via MALDIquant::intensityMatrix (as in the original examples)
#' ##IntMat <- MALDIquant::intensityMatrix(peaks)
#' ##rownames(IntMat) <- paste(CitrobacterRKImetadata$Strain_name_spot)
#' ##IntMat[is.na(IntMat)] <- 0
#' ##IntMat <- t(apply(IntMat, 1, function(x) x / max(x)))  # per-spectrum max norm
#' ##X <- t(IntMat)                                         # features in columns
#' ##Y <- factor(CitrobacterRKImetadata$Species)
#'
#' ## 3) Select discriminant m/z with "cvp" method
#' a <- MSclassifR::SelectionVar(
#'   X, Y,
#'   MethodSelection = "cvp",
#'   MethodValidation = "cv",
#'   PreProcessing = c("center","scale","nzv","corr"),
#'   NumberCV = 2,
#'   Metric = "Kappa"
#' )
#' sel_moz <- a$sel_moz
#'
#' ## 4) Train several models on the shortlisted m/z
#' model_lm  <- MSclassifR::LogReg(X = X, moz = sel_moz, Y = Y, number = 2,
#'  repeats = 2, Metric = "Kappa", kind = "linear")
#' model_nn  <- MSclassifR::LogReg(X = X, moz = sel_moz, Y = Y, number = 2,
#'  repeats = 2, Metric = "Kappa", kind = "nnet", Sampling = "up")
#' model_rf  <- MSclassifR::LogReg(X = X, moz = sel_moz, Y = Y, number = 2,
#'  repeats = 2, Metric = "Kappa", kind = "rf",  Sampling = "down")
#' model_xgb <- MSclassifR::LogReg(X = X, moz = sel_moz, Y = Y, number = 2,
#'  repeats = 2, Metric = "Kappa", kind = "xgb", Sampling = "smote")
#' model_svm <- MSclassifR::LogReg(X = X, moz = sel_moz, Y = Y, number = 2,
#'  repeats = 2, Metric = "Kappa", kind = "svm", Sampling = "up")
#'
#' Models <- list(
#'   model_lm$train_mod,
#'   model_nn$train_mod,
#'   model_rf$train_mod,
#'   model_xgb$train_mod,
#'   model_svm$train_mod
#' )
#'
#' ## 5) Predict classes for a subset of peaks; 6 Da tolerance for matching
#' prob_cat <- MSclassifR::PredictLogReg(
#'   peaks     = peaks[1:5],
#'   model     = Models,
#'   moz       = sel_moz,
#'   tolerance = 6,
#'   Reference = Y[1:5]
#' )
#' prob_cat
#'
#' ## 6) Meta-classifier strategy (several RF models + SMOTE + Fisher combine)
#' a2 <- MSclassifR::SelectionVar(X, Y, MethodSelection = "mda", Ntree = 5 * ncol(X))
#' sel_moz2 <- a2$sel_moz
#' models2 <- vector("list", 4L)
#' for (i in seq_along(models2)) {
#'   models2[[i]] <- MSclassifR::LogReg(
#'     X = X, moz = sel_moz2, Y = Y,
#'     number = 5, repeats = 5,
#'     kind = "rf", Metric = "Kappa",
#'     Sampling = "smote"
#'   )$train_mod
#' }
#' prob_cat2 <- MSclassifR::PredictLogReg(
#'   peaks = peaks,
#'   model = models2,
#'   moz   = sel_moz2,
#'   tolerance = 6,
#'   Reference = Y
#' )
#' }
#'
#' @references
#' Kuhn, M. (2008). Building predictive models in R using the caret package. Journal of Statistical Software, 28(1), 1–26.
#'
#' Alexandre Godmer, Yahia Benzerara, Emmanuelle Varon, Nicolas Veziris, Karen Druart,
#' Renaud Mozet, Mariette Matondo, Alexandra Aubry, Quentin Giai Gianetto (2025).
#' MSclassifR: An R package for supervised classification of mass spectra with machine learning methods.
#' Expert Systems with Applications, 294, 128796. doi:10.1016/j.eswa.2025.128796
#'
#' @seealso LogReg; SelectionVar; SelectionVarStat_fast; build_X_from_peaks_fast
#' @export
PredictLogReg <- function(peaks,
                            model,
                            moz,
                            tolerance = 6,
                            toleranceStep = 2,
                            normalizeFun = TRUE,
                            noMatch = 0,
                            Reference = NULL,
                            chunk_size = 10000L,
                            ncores = 1L,
                            verbose = FALSE) {

    # Accept a single model or list of models
    models <- if (!is.null(model$method)) list(model) else model
    moz <- sort(unique(as.numeric(moz)))

    # Build X from peaks
    X <- build_X_from_peaks_fast(
      peaks, moz, tolerance,
      normalize = normalizeFun,
      noMatch = noMatch,
      bump_if_empty = TRUE,
      toleranceStep = toleranceStep
    )

    # Spectrum names (as character)
    spec_name <- vapply(peaks, function(pk) {
      if (is.null(pk@metaData$fullName)) pk@metaData$file else pk@metaData$fullName
    }, character(1))
    spec_name <- as.character(spec_name)

    out_list <- vector("list", length(models))
    meth_names <- character(length(models))

    # Chunk rows for memory safety
    row_blocks <- split(seq_len(nrow(X)), ceiling(seq_len(nrow(X)) / chunk_size))

    # Helper: get class probabilities for a model for a data.frame block
    predict_probs <- function(mdl, newdata_df) {
      # Case A: caret::train object (supports type = "prob")
      if (inherits(mdl, "train")) {
        pr_try <- try(stats::predict(mdl, newdata_df, type = "prob"), silent = TRUE)
        if (!inherits(pr_try, "try-error")) {
          pr <- as.matrix(pr_try)
          colnames(pr) <- make.names(colnames(pr))
          pr[!is.finite(pr)] <- 0
          pr <- pmin(pmax(pr, 0), 1)
          return(pr)
        }
        # fallback: raw predictions -> one-hot
        raw <- stats::predict(mdl, newdata_df, type = "raw")
        levs <- if (!is.null(mdl$levels)) mdl$levels else levels(raw)
        levs <- make.names(levs)
        pr <- matrix(0, nrow = length(raw), ncol = length(levs), dimnames = list(NULL, levs))
        pr[cbind(seq_along(raw), match(make.names(as.character(raw)), levs))] <- 1
        return(pr)
      }
      # Case B: fast RF path — a list with $model = ranger object
      if (is.list(mdl) && !is.null(mdl$model) && inherits(mdl$model, "ranger")) {
        pr <- predict(mdl$model, data = newdata_df)$predictions
        pr <- as.matrix(pr)
        colnames(pr) <- make.names(colnames(pr))
        pr[!is.finite(pr)] <- 0
        pr <- pmin(pmax(pr, 0), 1)
        return(pr)
      }
      # Case C: a raw ranger model
      if (inherits(mdl, "ranger")) {
        pr <- predict(mdl, data = newdata_df)$predictions
        pr <- as.matrix(pr)
        colnames(pr) <- make.names(colnames(pr))
        pr[!is.finite(pr)] <- 0
        pr <- pmin(pmax(pr, 0), 1)
        return(pr)
      }
      # Fallback: try type = "prob"; else raw->one-hot
      pr_try <- try(stats::predict(mdl, newdata_df, type = "prob"), silent = TRUE)
      if (!inherits(pr_try, "try-error")) {
        pr <- as.matrix(pr_try)
        colnames(pr) <- make.names(colnames(pr))
        pr[!is.finite(pr)] <- 0
        pr <- pmin(pmax(pr, 0), 1)
        return(pr)
      }
      raw <- stats::predict(mdl, newdata_df, type = "raw")
      levs <- levels(factor(raw))
      levs <- make.names(levs)
      pr <- matrix(0, nrow = length(raw), ncol = length(levs), dimnames = list(NULL, levs))
      pr[cbind(seq_along(raw), match(make.names(as.character(raw)), levs))] <- 1
      pr
    }

    for (i in seq_along(models)) {
      mdl <- models[[i]]
      # Method name (caret::train has $method; fast RF path list has $method)
      if (!is.null(mdl$method)) {
        meth_names[i] <- paste0(mdl$method, "_", i)
      } else if (inherits(mdl, "ranger") || (is.list(mdl) && inherits(mdl$model, "ranger"))) {
        meth_names[i] <- paste0("ranger_", i)
      } else {
        meth_names[i] <- paste0("model_", i)
      }

      probs_accum <- NULL
      for (idx in row_blocks) {
        dfblock <- as.data.frame(X[idx, , drop = FALSE], check.names = FALSE)
        pr <- predict_probs(mdl, dfblock)

        # Align class columns across blocks within this model if needed
        if (is.null(probs_accum)) {
          probs_accum <- pr
        } else {
          if (!identical(colnames(probs_accum), colnames(pr))) {
            allc <- union(colnames(probs_accum), colnames(pr))
            expand <- function(M, allc) {
              out <- matrix(0, nrow = nrow(M), ncol = length(allc), dimnames = list(NULL, allc))
              out[, colnames(M)] <- M
              out
            }
            probs_accum <- expand(probs_accum, allc)
            pr <- expand(pr, allc)
          }
          probs_accum <- rbind(probs_accum, pr)
        }
      }
      out_list[[i]] <- data.frame(name = spec_name, method = meth_names[i], probs_accum, check.names = FALSE)
    }

    # Harmonize class columns across models before rbind
    class_cols_list <- lapply(out_list, function(df) setdiff(colnames(df), c("name", "method")))
    all_classes <- sort(Reduce(union, class_cols_list))

    out_list_aligned <- lapply(out_list, function(df) {
      missing <- setdiff(all_classes, colnames(df))
      for (m in missing) df[[m]] <- 0
      df <- df[, c("name", "method", all_classes), drop = FALSE]
      # force character columns
      df$name <- as.character(df$name)
      df$method <- as.character(df$method)
      df
    })

    Results <- do.call(rbind, out_list_aligned)
    rownames(Results) <- paste0(seq_len(nrow(Results)), ".")
    # ensure character types
    Results$name <- as.character(Results$name)
    Results$method <- as.character(Results$method)

    # Optional: multi-model combination
    if (length(models) > 1) {
      have_metap <- requireNamespace("metap", quietly = TRUE)
      cls_cols <- setdiff(colnames(Results), c("name", "method"))
      spl <- split(Results, Results$name)
      comb_rows <- lapply(names(spl), function(nm) {
        df <- spl[[nm]]
        probs <- as.matrix(df[, cls_cols, drop = FALSE])   # rows = models, cols = classes
        rows_to_bind <- list()

        # Fisher combine only if metap is available
        if (have_metap) {
          merge_p <- suppressWarnings(vapply(seq_len(ncol(probs)), function(j) {
            metap::sumlog(as.numeric(probs[, j]))$p
          }, numeric(1)))
          names(merge_p) <- colnames(probs)
          r1 <- c(name = nm, method = "comb_fisher", as.list(merge_p))
          rows_to_bind[[length(rows_to_bind) + 1L]] <- as.data.frame(t(r1), check.names = FALSE)
        }

        # Majority voting across models
        votes <- apply(probs, 1L, function(r) { v <- rep(0, ncol(probs)); v[which.max(r)] <- 1; v })
        if (is.vector(votes)) votes <- matrix(votes, nrow = ncol(probs), ncol = 1L)
        vote_frac <- rowMeans(votes)
        names(vote_frac) <- colnames(probs)
        r2 <- c(name = nm, method = "max_vote", as.list(vote_frac))
        rows_to_bind[[length(rows_to_bind) + 1L]] <- as.data.frame(t(r2), check.names = FALSE)

        do.call(rbind, rows_to_bind)
      })
      comb_df <- do.call(rbind, comb_rows)
      # ensure character types
      comb_df$name <- as.character(comb_df$name)
      comb_df$method <- as.character(comb_df$method)
      Results <- rbind(Results, comb_df)
    }

    # Predicted class per row
    class_cols <- setdiff(colnames(Results), c("name", "method"))
    Results$pred_max_p <- class_cols[
      max.col(as.matrix(Results[, class_cols, drop = FALSE]), ties.method = "first")
    ]

    # Confusion matrices if Reference provided (align by spectrum name)
    if (!is.null(Reference) && length(unique(as.character(Reference))) > 1) {
      # Named sanitized reference by spectrum name
      ref_vec <- setNames(make.names(as.character(Reference)), spec_name)
      levels_ref <- levels(factor(ref_vec))

      # Exclude combined rows
      df_eval <- Results[!(Results$method %in% c("comb_fisher", "max_vote")), , drop = FALSE]
      # ensure character
      df_eval$name <- as.character(df_eval$name)
      df_eval$method <- as.character(df_eval$method)

      df_methods <- split(df_eval, df_eval$method)

      Confusion.Matrix <- lapply(df_methods, function(dfm) {
        nm <- as.character(dfm$name)
        ref_ordered <- factor(ref_vec[nm], levels = levels_ref)
        preds <- factor(dfm$pred_max_p, levels = levels_ref)
        try(caret::confusionMatrix(preds, ref_ordered), silent = TRUE)
      })
      # keep only successful confusion matrices
      Confusion.Matrix <- Confusion.Matrix[!vapply(Confusion.Matrix, is.character, logical(1))]
      return(list(Prob.results = Results, Confusion.Matrix = Confusion.Matrix))
    }

    Results
  }

