
d_left_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1, distance_col = NULL) {
  # Validate method
  method <- match.arg(method, c("euclidean", "manhattan"))

  # Handle column specification
  if (is.null(by)) {
    by <- intersect(names(x), names(y))
  }
  if (length(by) == 0) {
    stop("No common column names found for joining")
  }

  # Create matrices for distance calculation
  x_matrix <- as.matrix(x[, by, drop = FALSE])
  y_matrix <- as.matrix(y[, by, drop = FALSE])

  # Calculate all pairwise distances using vectorized operations
  distances <- matrix(0, nrow = nrow(x), ncol = nrow(y))

  for (i in 1:length(by)) {
    # Extract vectors for current column
    x_col <- x_matrix[, i]
    y_col <- y_matrix[, i]

    # Calculate differences for all pairs
    if (method == "euclidean") {
      diff_matrix <- outer(x_col, y_col, function(a, b) (a - b)^2)
      distances <- distances + diff_matrix
    } else { # manhattan
      diff_matrix <- outer(x_col, y_col, function(a, b) abs(a - b))
      distances <- distances + diff_matrix
    }
  }

  # Finalize distance calculation
  if (method == "euclidean") {
    distances <- sqrt(distances)
  }

  # Find matches
  match_indices <- lapply(1:nrow(x), function(i) {
    matches <- which(distances[i, ] <= max_dist)
    if (length(matches) > 0) {
      # Take the match with minimum distance
      best_match <- matches[which.min(distances[i, matches])]
      return(data.frame(x_idx = i, y_idx = best_match,
                        distance = distances[i, best_match]))
    } else {
      return(NULL)
    }
  })

  # Combine match indices into a data frame
  match_df <- do.call(rbind, match_indices[!sapply(match_indices, is.null)])

  # Prepare the result dataframe with proper column handling

  # 1. Start with a copy of x dataframe
  result <- x

  # 2. Handle common columns (add .x suffix to original x columns)
  common_cols <- intersect(names(x), names(y))
  for (col in common_cols) {
    # Rename columns that appear in both dataframes
    # Avoid renaming the column if it's the distance column
    should_rename <- TRUE
    if (!is.null(distance_col)) {
      if (col == distance_col) {
        should_rename <- FALSE
      }
    }

    if (should_rename) {
      names(result)[names(result) == col] <- paste0(col, ".x")
    }
  }

  # 3. Create y columns
  if (nrow(match_df) > 0) {
    # Add all y columns (with .y suffix for common columns)
    for (col in names(y)) {
      if (col %in% common_cols) {
        # Common column gets .y suffix
        col_name <- paste0(col, ".y")
      } else {
        # Unique y column keeps its name
        col_name <- col
      }

      # Initialize with NA
      result[[col_name]] <- NA

      # Fill in matched values
      for (i in 1:nrow(match_df)) {
        x_idx <- match_df$x_idx[i]
        y_idx <- match_df$y_idx[i]
        result[x_idx, col_name] <- y[y_idx, col]
      }
    }
  } else {
    # No matches - still need to add y columns with NAs
    for (col in names(y)) {
      if (col %in% common_cols) {
        result[[paste0(col, ".y")]] <- NA
      } else {
        result[[col]] <- NA
      }
    }
  }

  # 4. Add distance column if requested
  if (!is.null(distance_col)) {
    result[[distance_col]] <- NA
    if (nrow(match_df) > 0) {
      for (i in 1:nrow(match_df)) {
        result[match_df$x_idx[i], distance_col] <- match_df$distance[i]
      }
    }
  }

  return(result)
}
