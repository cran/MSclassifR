#' Internal: choose a safe number of workers
#' @keywords internal
.safe_n_workers <- function(requested = NULL, fallback = 1L) {
  # default to fallback (1) unless a value is provided
  n <- if (is.null(requested)) as.integer(fallback) else as.integer(requested)
  if (is.na(n) || n < 1L) n <- 1L
  # Respect R CMD check core limit
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    n <- min(n, 2L)
  }
  n
}

`%||%` <- function(a, b) if (is.null(a)) b else a


