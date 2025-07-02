# Fast neighbor finding function with error checking
fast_find_neighbors <- function(data, nominal_indices, p_code, k) {
  n <- nrow(data)

  # Validate k is not larger than the dataset
  if (k >= n) {
    stop("Number of neighbors k must be less than the number of examples")
  }

  # Pre-compute distance matrix using efficient approach
  dist_matrix <- matrix(0, nrow = n, ncol = n)

  # Sequential version
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- calculate_distance(data[i,], data[j,], nominal_indices, p_code)
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    }
  }

  # Initialize neighbors matrix
  neighbors <- matrix(0, nrow = n, ncol = k)

  # Find k nearest neighbors for each example
  for (i in 1:n) {
    # Make a copy of distances and set self-distance to Inf
    distances <- dist_matrix[i, ]
    distances[i] <- Inf

    # Find k smallest distances
    nearest <- order(distances)[1:k]
    neighbors[i, ] <- nearest
  }

  return(neighbors)
}
