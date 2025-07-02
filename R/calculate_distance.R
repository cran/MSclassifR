# Optimized distance calculation function
calculate_distance <- function(x, y, nominal_indices, p_code) {
  # Safety check
  if (length(x) != length(y)) {
    stop("Vectors must have the same length for distance calculation")
  }

  # Fast distance calculation based on metric
  if (p_code >= 1) {
    # p-norm
    return(sum(abs(x - y)^p_code)^(1/p_code))
  } else if (p_code == 0) {
    # Chebyshev
    return(max(abs(x - y)))
  } else if (p_code == -1) {
    # Canberra
    sum_vals <- abs(x) + abs(y)
    # Avoid division by zero
    sum_vals[sum_vals == 0] <- 1
    return(sum(abs(x - y) / sum_vals))
  } else if (p_code == -2) {
    # Overlap (nominal only)
    if (length(nominal_indices) > 0) {
      return(sum(x[nominal_indices] != y[nominal_indices]))
    } else {
      return(0)
    }
  } else if (p_code == -3) {
    # HEOM - Heterogeneous Euclidean-Overlap Metric
    dist_num <- 0
    dist_nom <- 0

    # Calculate numeric distance
    if (length(nominal_indices) < length(x)) {
      numeric_indices <- setdiff(1:length(x), nominal_indices)
      dist_num <- sum((x[numeric_indices] - y[numeric_indices])^2)
    }

    # Calculate nominal distance
    if (length(nominal_indices) > 0) {
      dist_nom <- sum(x[nominal_indices] != y[nominal_indices])
    }

    return(sqrt(dist_num + dist_nom))
  } else if (p_code == -4) {
    # HVDM - simplified version
    dist_num <- 0
    dist_nom <- 0

    # Calculate numeric distance (simplified)
    if (length(nominal_indices) < length(x)) {
      numeric_indices <- setdiff(1:length(x), nominal_indices)
      dist_num <- sum((x[numeric_indices] - y[numeric_indices])^2)
    }

    # Calculate nominal distance (simplified)
    if (length(nominal_indices) > 0) {
      dist_nom <- sum(x[nominal_indices] != y[nominal_indices])
    }

    return(sqrt(dist_num + dist_nom))
  }

  # Default Euclidean distance
  return(sqrt(sum((x - y)^2)))
}
