#' Parallel helpers (internal)
#' @keywords internal
#' @noRd
.limited_cores <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) n <- 1L
  # CRAN sets this env var; respect 2-core cap when TRUE
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  nmax <- if (identical(lim, "TRUE")) {
    2L
  } else {
    cores <- parallel::detectCores()
    if (is.na(cores)) 1L else cores
  }
  max(1L, min(n, nmax))
}

#' Safe worker count (internal)
#' @keywords internal
#' @noRd
.safe_n_workers <- function(n_workers, fallback = 1L) {
  n_req <- if (is.null(n_workers) || !is.finite(n_workers) || n_workers < 1L) {
    fallback
  } else {
    as.integer(n_workers)
  }
  .limited_cores(n_req)
}
