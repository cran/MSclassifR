SelectionVar <- function(X, Y,
                         MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA", "mda", "cvp", "boruta"),
                         MethodValidation = c("cv", "repeatedcv", "LOOCV"),
                         PreProcessing = c("center", "scale", "nzv", "corr"),
                         Metric = c("Kappa", "Accuracy"),
                         Sampling = c("no", "up","down", "smote"),
                         NumberCV = NULL, RepeatsCV = NULL,
                         Sizes,
                         Ntree = 1000,
                         ncores = 2,
                         threshold = 0.01,
                         ncomp.max = 10,
                         nbf = 0) {
  # Silence codetools NOTE for sPLSDA branch
  cp <- NULL
  UniqueVariables_keep <- NULL

  # Cap cores for CRAN and safety
  ncores_eff <- .limited_cores(ncores)

  # 0) Make inputs modeling-safe
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (anyNA(X)) X[is.na(X)] <- 0
  Y <- factor(Y)

  # 1) Match user choices
  method <- match.arg(MethodSelection)
  methodValid <- match.arg(MethodValidation)
  Metric <- match.arg(Metric)
  SamplingM <- match.arg(Sampling)

  # 2) Sensible defaults for CV
  if (methodValid == "cv") {
    if (is.null(NumberCV)) NumberCV <- 5L
    RepeatsCV <- NULL
  } else if (methodValid == "repeatedcv") {
    if (is.null(NumberCV)) NumberCV <- 5L
    if (is.null(RepeatsCV)) RepeatsCV <- 3L
  } else { # LOOCV
    NumberCV <- NULL
    RepeatsCV <- NULL
  }

  # 3) Sampling
  switch(SamplingM,
         "no" = { message("No sampling method selected") },
         "up" = {
           message("Up sampling method selected")
           upTrain <- caret::upSample(x = X, y = Y, list = TRUE)
           X <- upTrain$x; Y <- factor(upTrain$y)
         },
         "down" = {
           message("Down sampling method selected")
           downTrain <- caret::downSample(x = X, y = Y, list = TRUE)
           X <- downTrain$x; Y <- factor(downTrain$y)
         },
         "smote" = {
           message("Smote sampling method selected (internal smote_classif)")
           mozv <- colnames(X)
           dataSMOTE <- data.frame(Y = Y, X, check.names = FALSE)
           Smoted <- smote_classif(Y ~ ., dataSMOTE, C.perc = "balance")
           X <- as.matrix(Smoted[, -1, drop = FALSE]); X[is.na(X)] <- 0; colnames(X) <- mozv
           Y <- factor(Smoted$Y)
         }
  )

  # 4) Methods
  results <- switch(method,

                    "RFERF" = {
                      if (missing(Sizes) || length(Sizes) == 0L)
                        stop("Sizes must be provided for RFERF.", call. = FALSE)
                      message("Selection variables with RFE and RF method")
                      control <- caret::rfeControl(functions = caret::rfFuncs,
                                                   method = methodValid,
                                                   number = NumberCV,
                                                   repeats = RepeatsCV)
                      resultmodel <- caret::rfe(X, Y,
                                                rfeControl = control,
                                                sizes = Sizes,
                                                preProc = PreProcessing,
                                                metric = Metric)
                      nbv <- if (Metric == "Accuracy") {
                        resultmodel$results$Variables[which.max(resultmodel$results$Accuracy[1:(nrow(resultmodel$results) - 1)])]
                      } else {
                        resultmodel$results$Variables[which.max(resultmodel$results$Kappa[1:(nrow(resultmodel$results) - 1)])]
                      }
                      smoz <- unique(resultmodel$variables$var[resultmodel$variables$Variables == nbv])
                      list(result = resultmodel, sel_moz = sort(smoz))
                    },

                    "RFEGlmnet" = {
                      if (missing(Sizes) || length(Sizes) == 0L)
                        stop("Sizes must be provided for RFEGlmnet.", call. = FALSE)
                      message("Selection variables with RFE and Glmnet method")
                      control <- caret::rfeControl(functions = caret::caretFuncs,
                                                   method = methodValid,
                                                   number = NumberCV,
                                                   repeats = RepeatsCV)
                      MyPerf.rfeglmnet <- try(
                        resultmodel <- caret::rfe(X, Y,
                                                  rfeControl = control,
                                                  sizes = Sizes,
                                                  method = "glmnet",
                                                  preProc = PreProcessing,
                                                  metric = Metric),
                        silent = TRUE
                      )
                      if (inherits(MyPerf.rfeglmnet, "try-error")) {
                        warning("glmnet: fallback to LOOCV due to an error.")
                        control <- caret::rfeControl(functions = caret::caretFuncs, method = "LOOCV")
                        resultmodel <- caret::rfe(X, Y,
                                                  rfeControl = control,
                                                  sizes = Sizes,
                                                  method = "glmnet",
                                                  preProc = PreProcessing,
                                                  metric = Metric)
                      }
                      nbv <- if (Metric == "Accuracy") {
                        resultmodel$results$Variables[which.max(resultmodel$results$Accuracy[1:(nrow(resultmodel$results) - 1)])]
                      } else {
                        resultmodel$results$Variables[which.max(resultmodel$results$Kappa[1:(nrow(resultmodel$results) - 1)])]
                      }
                      smoz <- unique(resultmodel$variables$var[resultmodel$variables$Variables == nbv])
                      list(result = resultmodel, sel_moz = sort(smoz))
                    },

                    "sPLSDA" = {
                      if (missing(Sizes) || length(Sizes) == 0L)
                        stop("Sizes must be provided (test.keepX grid) for sPLSDA.", call. = FALSE)
                      message("Selection variables with sPLSDA method")

                      # Map validation for mixOmics
                      vali <- if (methodValid == "LOOCV") "loo" else "Mfold"
                      nrep <- if (methodValid == "repeatedcv") if (is.null(RepeatsCV)) 1L else RepeatsCV else 1L

                      SCALEm <- "scale" %in% PreProcessing
                      NZVm <- "nzv" %in% PreProcessing
                      if ("corr" %in% PreProcessing) warning("The 'corr' preprocessing method is not used with sPLSDA.")
                      if ("center" %in% PreProcessing) warning("The 'center' preprocessing method is not used with sPLSDA.")

                      ncomp.max <- min(ncomp.max, max(1L, nrow(X) - 1L))
                      MyResult.plsda <- mixOmics::plsda(X, Y, ncomp = ncomp.max, scale = SCALEm, near.zero.var = NZVm)

                      MyPerf.plsda <- try(
                        mixOmics::perf(MyResult.plsda, validation = vali, folds = NumberCV, nrepeat = nrep, progressBar = FALSE),
                        silent = TRUE
                      )
                      if (inherits(MyPerf.plsda, "try-error")) {
                        warning("plsda: fallback to LOOCV due to singular system.")
                        ok <- FALSE
                        while (!ok && ncomp.max > 1) {
                          ncomp.max <- ncomp.max - 1
                          MyResult.plsda_k <- mixOmics::plsda(X, Y, ncomp = ncomp.max, scale = SCALEm, near.zero.var = NZVm)
                          tst <- try(mixOmics::perf(MyResult.plsda_k, validation = "loo",
                                                    folds = NumberCV, nrepeat = 1L, progressBar = FALSE),
                                     silent = TRUE)
                          ok <- !inherits(tst, "try-error")
                          if (ok) MyPerf.plsda <- tst
                        }
                      }

                      adiff <- abs(diff(as.numeric(MyPerf.plsda$error.rate$BER[, 1]))) > threshold
                      k <- 1L
                      while (k <= length(adiff) && isTRUE(adiff[k])) k <- k + 1L
                      ncomp <- max(2L, k)

                      tune.splsda.srbct <- try(
                        mixOmics::tune.splsda(X, Y, ncomp = ncomp, validation = vali, folds = NumberCV,
                                              dist = "max.dist", progressBar = FALSE, measure = "BER",
                                              test.keepX = Sizes),
                        silent = TRUE
                      )
                      if (inherits(tune.splsda.srbct, "try-error")) {
                        warning("splsda: fallback to LOOCV during tuning due to singular system.")
                        tune.splsda.srbct <- mixOmics::tune.splsda(X, Y, ncomp = ncomp, validation = "loo", folds = NumberCV,
                                                                   dist = "max.dist", progressBar = FALSE, measure = "BER",
                                                                   test.keepX = Sizes)
                      }

                      select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]
                      splsda.train <- mixOmics::splsda(X, Y, scale = SCALEm, ncomp = ncomp, keepX = select.keepX, near.zero.var = NZVm)

                      compo <- seq_len(ncomp)
                      cp <- vector("list", length(compo))
                      variables_keep <- vector("list", length(compo))
                      for (i in seq_along(compo)) {
                        cp[[i]] <- mixOmics::plotLoadings(splsda.train, comp = compo[[i]],
                                                          title = paste("comp", compo[[i]], sep = "_"),
                                                          contrib = "max", method = "mean", plot = FALSE)
                        variables_keep[[i]] <- data.frame(
                          X_var = row.names(cp[[i]]),
                          Group = cp[[i]]$GroupContrib,
                          Importance = cp[[i]]$importance,
                          stringsAsFactors = FALSE
                        )
                      }
                      names(cp) <- paste0("comp_raw_", compo)
                      names(variables_keep) <- paste0("comp_", compo)
                      UniqueVariables_keep <- do.call("rbind", variables_keep)
                      UniqueVariables_keep <- UniqueVariables_keep[!duplicated(UniqueVariables_keep[, 1]), , drop = FALSE]

                      list(Raw_data = cp,
                           selected_variables = UniqueVariables_keep,
                           sel_moz = sort(UniqueVariables_keep$X_var))
                    },

                    "VSURF" = {
                      message("Selection variables with VSURF method")
                      resultsmodel <- VSURF::VSURF(x = X, y = Y, ntree = Ntree,
                                                   nfor.thres = 50, nfor.interp = 50, nfor.pred = 50, nsd = 100)
                      sel_moz <- colnames(X)[resultsmodel[["varselect.pred"]]]
                      list(result = resultsmodel, sel_moz = sort(sel_moz))
                    },

                    "mda" = {
                      message("Selection variables with mda method")
                      if (requireNamespace("ranger", quietly = TRUE)) {
                        return(fast_mda(X, Y, ntree = Ntree, nbf = nbf, nthreads = ncores_eff, seed = 123))
                      }
                      if (!requireNamespace("randomForest", quietly = TRUE))
                        stop("Method 'mda' requires 'ranger' (preferred) or 'randomForest'.", call. = FALSE)

                      # Fallback: randomForest (single-threaded)
                      if (nbf > 0) {
                        X0 <- matrix(stats::runif(nbf * nrow(X), min = min(X), max = max(X)), nrow(X), nbf)
                        colnames(X0) <- paste0("false_", seq_len(ncol(X0)))
                        Xn <- cbind(X, X0)
                      } else Xn <- X

                      rf <- randomForest::randomForest(x = Xn, y = Y, ntree = Ntree,
                                                       mtry = max(1L, floor(sqrt(ncol(X)))), importance = TRUE, keep.forest = FALSE)
                      vi_mat <- randomForest::importance(rf, type = 1, scale = FALSE)
                      vi_vec <- if (is.matrix(vi_mat)) setNames(as.numeric(vi_mat[, ncol(vi_mat)]), rownames(vi_mat)) else vi_mat
                      vi_vec[!is.finite(vi_vec)] <- 0

                      if (nbf > 0) {
                        vi_true <- vi_vec[colnames(X)]
                        vi_false_neg <- vi_vec[!names(vi_vec) %in% colnames(X)]
                        vi_false_neg <- vi_false_neg[is.finite(vi_false_neg) & vi_false_neg < 0]
                        vi1 <- c(vi_true, vi_false_neg)
                      } else {
                        vi1 <- vi_vec[colnames(X)]
                      }

                      imp_neg <- vi1[vi1 < 0]
                      if (length(imp_neg) == 0L) {
                        pi0f <- 0
                      } else {
                        imp_null <- c(imp_neg, -imp_neg)
                        q_ext <- seq(0.75, 1, by = 0.01)
                        Fall <- stats::ecdf(vi1)
                        pi0_raw <- vapply(q_ext, function(q) {
                          qin <- stats::quantile(imp_null, q, na.rm = TRUE); min(Fall(qin) / q, 1)
                        }, numeric(1))
                        if (nbf > 0) {
                          Nfn <- sum(vi_false_neg < 0)
                          pi0f <- (min(pi0_raw) * (ncol(X) + Nfn) - Nfn) / ncol(X)
                        } else {
                          pi0f <- min(pi0_raw)
                        }
                      }

                      nb_to_sel <- max(1L, floor(ncol(X) * (1 - pi0f)))
                      vi_true_only <- vi1[colnames(X)]
                      sel_moz <- names(vi_true_only)[order(-vi_true_only)][seq_len(nb_to_sel)]
                      imp_sel <- vi_true_only[sel_moz]
                      list(nb_to_sel = nb_to_sel, sel_moz = sel_moz, imp_sel = imp_sel)
                    },

                    "cvp" = {
                      message("Selection variables with cvp method")
                      if (is.null(NumberCV)) NumberCV <- 5L
                      if (requireNamespace("ranger", quietly = TRUE)) {
                        return(fast_cvpvi(X, Y, k = NumberCV, ntree = Ntree, nbf = nbf, nthreads = ncores_eff, seed = 123))
                      }
                      if (!requireNamespace("vita", quietly = TRUE))
                        stop("Method 'cvp' requires 'ranger' (preferred) or 'vita'.", call. = FALSE)

                      # Fallback: vita::CVPVI (cap ncores)
                      if (nbf > 0) {
                        X0 <- matrix(stats::runif(nbf * nrow(X), min = min(X), max = max(X)), nrow(X), nbf)
                        colnames(X0) <- paste0("false_", seq_len(ncol(X0)))
                        Xn <- cbind(X, X0)
                      } else Xn <- X

                      cv_vi <- vita::CVPVI(Xn, as.numeric(Y), k = NumberCV, ntree = Ntree, ncores = ncores_eff)
                      vi <- as.numeric(cv_vi$cv_varim[, 1L]); names(vi) <- rownames(cv_vi$cv_varim)
                      vi[!is.finite(vi)] <- 0

                      if (nbf > 0) {
                        vi_true <- vi[colnames(X)]
                        vi_false_neg <- vi[!names(vi) %in% colnames(X)]
                        vi_false_neg <- vi_false_neg[is.finite(vi_false_neg) & vi_false_neg < 0]
                        vi1 <- c(vi_true, vi_false_neg)
                      } else {
                        vi1 <- vi[colnames(X)]
                      }

                      imp_neg <- vi1[vi1 < 0]
                      if (length(imp_neg) == 0L) {
                        pi0f <- 0
                      } else {
                        imp_null <- c(imp_neg, -imp_neg)
                        q_ext <- seq(0.75, 1, by = 0.01)
                        Fall <- stats::ecdf(vi1)
                        pi0_raw <- vapply(q_ext, function(q) {
                          qin <- stats::quantile(imp_null, q, na.rm = TRUE); min(Fall(qin) / q, 1)
                        }, numeric(1))
                        if (nbf > 0) {
                          Nfn <- sum(vi_false_neg < 0)
                          pi0f <- (min(pi0_raw) * (ncol(X) + Nfn) - Nfn) / ncol(X)
                        } else {
                          pi0f <- min(pi0_raw)
                        }
                      }

                      nb_to_sel <- max(1L, floor(ncol(X) * (1 - pi0f)))
                      vi_true_only <- vi1[colnames(X)]
                      sel_moz <- names(vi_true_only)[order(-vi_true_only)][seq_len(nb_to_sel)]
                      imp_sel <- vi_true_only[sel_moz]
                      list(nb_to_sel = nb_to_sel, sel_moz = sel_moz, imp_sel = imp_sel)
                    },

                    "boruta" = {
                      message("Selection variables with Boruta method")
                      e <- Boruta::Boruta(x = X, y = Y, maxRuns = 3 * ncol(X))
                      sel_moz <- colnames(X)[e$finalDecision == "Confirmed"]
                      list(nb_to_sel = length(sel_moz), sel_moz = sel_moz)
                    }
  )

  return(results)
}
