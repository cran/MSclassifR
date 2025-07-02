# Optimized synthetic example generation function
fast_generate_synthetic <- function(dat, k, n, p_code) {
  # Extract dimensions
  n_examples <- nrow(dat)
  n_features <- ncol(dat)
  target_col <- n_features

  # Handle special case of no examples to generate
  if (n <= 0) return(dat[0, ])

  # Separate numeric and factor columns for efficiency
  is_factor <- sapply(dat, is.factor)
  # Keep track of nominal attribute indices
  nomatr <- which(is_factor[-target_col])

  # Create numeric matrix from data (much faster for distance calculation)
  T <- matrix(0, nrow = n_examples, ncol = n_features - 1)
  for (j in 1:(n_features - 1)) {
    if (is_factor[j]) {
      T[, j] <- as.integer(dat[[j]])
    } else {
      T[, j] <- as.numeric(dat[[j]])
    }
  }

  # Fast distance calculation using optimized functions
  kNNs <- fast_find_neighbors(T, nomatr, p_code, k)

  # Safety check for kNNs dimensions
  if (nrow(kNNs) != n_examples || ncol(kNNs) != k) {
    stop(sprintf("Error in nearest neighbors calculation. Expected dimensions [%d, %d], got [%d, %d]",
                 n_examples, k, nrow(kNNs), ncol(kNNs)))
  }

  # Calculate synthetic examples to generate per example
  nexs <- floor(n / n_examples)
  extra <- n - nexs * n_examples

  # Pre-allocate result matrix (much faster than growing)
  newM <- matrix(0, nrow = n, ncol = n_features - 1)

  # Generate bulk of synthetic examples
  if (nexs > 0) {
    # Vector of random neighbor selections (faster than per-iteration sampling)
    random_neighbors <- sample(1:k, nexs * n_examples, replace = TRUE)
    random_weights <- runif(nexs * n_examples)

    row_idx <- 0
    for (i in 1:n_examples) {
      for (j in 1:nexs) {
        row_idx <- row_idx + 1
        neig_idx <- random_neighbors[row_idx]
        weight <- random_weights[row_idx]

        # Get neighbor index with safety check
        neig <- kNNs[i, neig_idx]

        # Fast difference calculation
        difs <- T[neig, ] - T[i, ]
        newM[row_idx, ] <- T[i, ] + weight * difs

        # Handle nominal attributes
        if (length(nomatr) > 0) {
          for (a in nomatr) {
            # Random selection between values
            if (runif(1) > 0.5) {
              newM[row_idx, a] <- T[neig, a]
            } else {
              newM[row_idx, a] <- T[i, a]
            }
          }
        }
      }
    }
  }

  # Generate extra examples if needed
  if (extra > 0) {
    idx <- sample(1:n_examples, extra)
    random_neighbors <- sample(1:k, extra, replace = TRUE)
    random_weights <- runif(extra)

    for (i in 1:extra) {
      example_idx <- idx[i]
      neig_idx <- random_neighbors[i]
      weight <- random_weights[i]

      # Get neighbor index with safety check
      neig <- kNNs[example_idx, neig_idx]

      # Fast difference calculation
      difs <- T[neig, ] - T[example_idx, ]
      newM[nexs * n_examples + i, ] <- T[example_idx, ] + weight * difs

      # Handle nominal attributes
      if (length(nomatr) > 0) {
        for (a in nomatr) {
          # Random selection between values
          if (runif(1) > 0.5) {
            newM[nexs * n_examples + i, a] <- T[neig, a]
          } else {
            newM[nexs * n_examples + i, a] <- T[example_idx, a]
          }
        }
      }
    }
  }

  # Convert results to data frame efficiently
  newCases <- as.data.frame(newM)

  # Set appropriate column names
  colnames(newCases) <- names(dat)[-target_col]

  # Restore factor levels for nominal attributes efficiently
  for (j in which(is_factor[-target_col])) {
    # Make sure we're converting valid factor levels
    factor_values <- newM[, j]
    valid_levels <- sort(unique(factor_values))

    newCases[[j]] <- factor(factor_values,
                            levels = valid_levels,
                            labels = levels(dat[[j]])[valid_levels])
  }

  # Add target column
  newCases[[target_col]] <- factor(rep(dat[[target_col]][1], nrow(newCases)),
                                   levels = levels(dat[[target_col]]))
  names(newCases)[target_col] <- names(dat)[target_col]

  return(newCases)
}
