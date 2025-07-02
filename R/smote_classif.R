smote_classif <- function(form, dat, C.perc = "balance", k = 5, repl = FALSE,
                          dist = "Euclidean", p = 2) {
  # Fast SMOTE implementation for classification problems
  #
  # Args:
  #   form: a model formula
  #   dat: the original training set (with the unbalanced distribution)
  #   C.perc: named list with percentages of under/over-sampling, or "balance"/"extreme"
  #   k: number of neighbors to consider
  #   repl: logical, whether to perform sampling with replacement
  #   dist: distance measure ("Euclidean", "Manhattan", "Chebyshev", etc.)
  #   p: parameter used when a p-norm is computed
  #
  # Returns: a new data frame modified through the SMOTE algorithm

  # Input validation
  if (any(is.na(dat))) {
    stop("The data set provided contains NA values!")
  }

  # Find target variable position
  tgt <- which(names(dat) == as.character(form[[2]]))
  if (length(tgt) == 0) {
    stop("Target variable not found in data")
  }

  # Extract target and convert to factor if needed
  y_col <- dat[[tgt]]
  if (!is.factor(y_col)) y_col <- as.factor(y_col)
  dat[[tgt]] <- y_col

  # Get class information
  class_names <- levels(y_col)
  class_counts <- tabulate(as.integer(y_col))
  names(class_counts) <- class_names

  # Handle columns rearrangement efficiently
  needs_reordering <- (tgt < ncol(dat))
  if (needs_reordering) {
    original_order <- names(dat)
    # Move target to end (more efficient for processing)
    dat <- dat[, c(setdiff(seq_len(ncol(dat)), tgt), tgt)]
  }

  # Determine sampling strategy (vectorized operations)
  if (is.list(C.perc)) {
    # User-specified sampling
    sampling_factors <- rep(1, length(class_names))
    names(sampling_factors) <- class_names

    for (cls in names(C.perc)) {
      if (cls %in% class_names) sampling_factors[cls] <- C.perc[[cls]]
    }
  } else if (C.perc == "balance") {
    # Balance classes
    target_size <- round(sum(class_counts)/length(class_counts), 0)
    sampling_factors <- target_size / class_counts
    names(sampling_factors) <- class_names
  } else if (C.perc == "extreme") {
    # Extreme balancing
    med <- sum(class_counts)/length(class_counts)
    target_sizes <- round(med^2/class_counts * sum(class_counts)/sum(med^2/class_counts), 0)
    sampling_factors <- target_sizes / class_counts
    names(sampling_factors) <- class_names
  } else {
    stop("Please provide a list with classes to under-/over-sample or indicate 'balance' or 'extreme'.")
  }

  # Use list to collect results (more efficient than rbind in a loop)
  result_chunks <- list()
  chunk_index <- 1

  # Fast categorization of classes
  under_classes <- names(sampling_factors)[sampling_factors < 1]
  over_classes <- names(sampling_factors)[sampling_factors > 1]
  same_classes <- names(sampling_factors)[sampling_factors == 1]

  # Process unchanged classes (vectorized)
  if (length(same_classes) > 0) {
    same_indices <- which(dat[[ncol(dat)]] %in% same_classes)
    if (length(same_indices) > 0) {
      result_chunks[[chunk_index]] <- dat[same_indices, , drop = FALSE]
      chunk_index <- chunk_index + 1
    }
  }

  # Process undersampled classes (vectorized where possible)
  if (length(under_classes) > 0) {
    for (cls in under_classes) {
      cls_indices <- which(dat[[ncol(dat)]] == cls)
      if (length(cls_indices) > 0) {
        sample_size <- max(1, round(sampling_factors[cls] * length(cls_indices)))
        selected <- sample(cls_indices, sample_size, replace = repl)
        result_chunks[[chunk_index]] <- dat[selected, , drop = FALSE]
        chunk_index <- chunk_index + 1
      }
    }
  }

  # Map distance metric to code (more efficient switch)
  p_code <- switch(dist,
                   "Chebyshev" = 0,
                   "Manhattan" = 1,
                   "Euclidean" = 2,
                   "Canberra" = -1,
                   "Overlap" = -2,
                   "HEOM" = -3,
                   "HVDM" = -4,
                   "p-norm" = p,
                   stop("Distance measure not available!"))

  # Process oversampled classes (optimized)
  if (length(over_classes) > 0) {
    for (cls in over_classes) {
      cls_indices <- which(dat[[ncol(dat)]] == cls)
      n_minority <- length(cls_indices)

      if (n_minority == 0) {
        next
      } else if (n_minority == 1) {
        # One example case - create replicas (vectorized)
        warning(paste("SmoteClassif: Unable to use SMOTE with 1 example of class", cls,
                      "- creating replicas instead."), call. = FALSE)

        n_replicas <- ceiling(sampling_factors[cls] - 1)
        if (n_replicas > 0) {
          replicated_data <- dat[rep(cls_indices, n_replicas), , drop = FALSE]
          result_chunks[[chunk_index]] <- replicated_data
          chunk_index <- chunk_index + 1
        }

        # Add original
        result_chunks[[chunk_index]] <- dat[cls_indices, , drop = FALSE]
        chunk_index <- chunk_index + 1

      } else if (n_minority <= k) {
        # Small minority class case
        adjusted_k <- n_minority - 1
        warning(paste("SmoteClassif: Class", cls, "has only", n_minority,
                      "examples. Using k =", adjusted_k, "for nearest neighbors."), call. = FALSE)

        minority_data <- dat[cls_indices, , drop = FALSE]

        # Generate synthetic examples
        to_generate <- ceiling((sampling_factors[cls] - 1) * n_minority)
        if (to_generate > 0) {
          synthetic_data <- fast_generate_synthetic(
            minority_data,
            adjusted_k,
            to_generate,
            p_code
          )

          result_chunks[[chunk_index]] <- synthetic_data
          chunk_index <- chunk_index + 1
        }

        # Add original examples
        result_chunks[[chunk_index]] <- minority_data
        chunk_index <- chunk_index + 1

      } else {
        # Normal case - enough examples
        minority_data <- dat[cls_indices, , drop = FALSE]

        # Generate synthetic examples
        to_generate <- ceiling((sampling_factors[cls] - 1) * n_minority)
        if (to_generate > 0) {
          synthetic_data <- fast_generate_synthetic(
            minority_data,
            k,
            to_generate,
            p_code
          )

          result_chunks[[chunk_index]] <- synthetic_data
          chunk_index <- chunk_index + 1
        }

        # Add original examples
        result_chunks[[chunk_index]] <- minority_data
        chunk_index <- chunk_index + 1
      }
    }
  }

  # Combine results efficiently
  if (length(result_chunks) == 0) {
    warning("No data generated. Check your parameters.")
    return(dat)
  }

  newdata <- do.call(rbind, result_chunks)

  # Restore original column order if needed
  if (needs_reordering) {
    newdata <- newdata[, original_order]
  }

  return(newdata)
}
