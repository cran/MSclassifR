#' SMOTE for classification datasets
#'
#' Generate synthetic examples for minority classes using the SMOTE idea,
#' to balance a classification dataset.
#'
#' The function supports multi-class data. With strategy = "balance" (default),
#' each class is oversampled up to the size of the largest class. With
#' strategy = "perc", each class c is oversampled by round(n_c * perc/100).
#' Neighbors are computed within each class.
#'
#' @param formula A model formula target ~ predictors indicating the response and predictors.
#' @param data A data.frame containing the variables in the model.
#' @param k Integer, number of nearest neighbors used by SMOTE (default 5).
#' @param strategy One of "balance" (oversample to the max class size) or "perc"
#'   (oversample each class by a percentage). Default "balance".
#' @param perc Numeric percentage used when strategy = "perc" (e.g., 100 means
#'   generate as many synthetic examples as existing in the class). Ignored for "balance".
#' @param metric Distance metric for neighbor search: one of
#'   "euclidean", "manhattan", "chebyshev", "canberra", "overlap", "heom", "hvdm", "pnorm".
#'   Default "euclidean".
#' @param p Numeric p for the p-norm when metric = "pnorm"; also used implicitly for
#'   "euclidean" (p=2) and "manhattan" (p=1). Default 2.
#' @param seed Optional integer seed for reproducibility.
#' @param C.perc Deprecated. Backward-compatibility alias for oversampling
#'   control. If character "balance", mapped to strategy = "balance".
#'   If a single numeric, mapped to strategy = "perc" and perc = C.perc.
#'   Other forms are ignored with a warning.
#'
#' @return A data.frame with synthetic rows appended, same columns and types as input.
#'
#' @examples
#' \donttest{
#' data(iris)
#' imbal_iris <- iris[c(1:40, 51:100, 101:110), ]
#' table(imbal_iris$Species)
#' balanced_iris <- smote_classif(Species ~ ., imbal_iris)
#' table(balanced_iris$Species)
#' }
#'
#' @export
smote_classif <- function(formula, data,
                          k = 5,
                          strategy = c("balance", "perc"),
                          perc = NULL,
                          metric = c("euclidean", "manhattan", "chebyshev", "canberra",
                                     "overlap", "heom", "hvdm", "pnorm"),
                          p = 2,
                          seed = NULL,
                          C.perc = NULL) {
  # Map deprecated C.perc to new args if present
  if (!is.null(C.perc)) {
    if (is.character(C.perc) && length(C.perc) == 1 && identical(tolower(C.perc), "balance")) {
      strategy <- "balance"
    } else if (is.numeric(C.perc) && length(C.perc) == 1) {
      strategy <- "perc"
      perc <- as.numeric(C.perc)
    } else {
      warning("C.perc is deprecated; use strategy='balance' or strategy='perc' + perc=<number>.")
    }
  }

  strategy <- match.arg(strategy)
  metric   <- match.arg(metric)
  if (!is.null(seed)) set.seed(seed)

  # Build model frame (drops rows with NA)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  if (ncol(mf) < 2) stop("Formula must be of the form target ~ predictors.")

  # Target and predictors
  target_name <- names(mf)[1]
  y <- mf[[1]]
  if (!is.factor(y)) y <- factor(y)
  X <- mf[-1]

  # Coerce character predictors to factors
  is_char <- vapply(X, is.character, logical(1))
  if (any(is_char)) X[is_char] <- lapply(X[is_char], factor)

  # Metric mapping to p_code
  p_code <- switch(metric,
                   euclidean = 2,
                   manhattan = 1,
                   chebyshev = 0,
                   canberra  = -1,
                   overlap   = -2,
                   heom      = -3,
                   hvdm      = -4,
                   pnorm     = if (p <= 0) stop("p must be > 0 for pnorm") else p)

  # Class sizes
  tbl <- table(y)
  max_n <- max(tbl)

  # Helper to generate synthetic rows for one class
  synth_one_class <- function(cls, n_new) {
    if (n_new <= 0) return(mf[0, , drop = FALSE])

    idx <- which(y == cls)
    df_c <- X[idx, , drop = FALSE]
    y_c  <- y[idx, drop = TRUE]
    n_c <- length(idx)

    if (n_c <= 1) {
      samp <- df_c[rep(1, n_new), , drop = FALSE]
      out <- cbind(samp, setNames(data.frame(factor(rep(cls, n_new), levels = levels(y))), target_name))
      return(out)
    }

    k_use <- min(k, n_c - 1)
    if (k_use < 1) {
      samp <- df_c[sample.int(n_c, n_new, replace = TRUE), , drop = FALSE]
      out <- cbind(samp, setNames(data.frame(factor(rep(cls, n_new), levels = levels(y))), target_name))
      return(out)
    }

    dat_c <- df_c
    dat_c[[target_name]] <- y_c

    syn <- fast_generate_synthetic(dat = dat_c, k = k_use, n = n_new, p_code = p_code)

    if (any(is_char)) {
      for (nm in names(X)[is_char]) syn[[nm]] <- as.character(syn[[nm]])
    }

    syn <- syn[c(names(X), target_name)]
    syn
  }

  # How many to generate per class
  add_list <- vector("list", length = length(levels(y)))
  names(add_list) <- levels(y)

  for (cls in levels(y)) {
    n_c <- as.integer(tbl[cls])
    n_new <- if (strategy == "balance") {
      max_n - n_c
    } else {
      if (is.null(perc)) stop("Provide 'perc' when strategy = 'perc'.")
      round(n_c * perc / 100)
    }
    add_list[[cls]] <- synth_one_class(cls, n_new)
  }

  synthetic <- do.call(rbind, add_list)
  out <- rbind(cbind(X, setNames(data.frame(y), target_name)), synthetic)
  rownames(out) <- NULL
  out
}
