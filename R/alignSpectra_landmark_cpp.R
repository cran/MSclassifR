#' Internal: fast landmark-based alignment via Rcpp (wrapper)
#'
#' Applies cpp_match_anchors and cpp_warp_mz_piecewise per spectrum to warp
#' the m/z axis to a reference using piecewise-linear mapping through matched anchors.
#' Not exported; used by SignalProcessingUltra when alignSpectra_method = "landmark_cpp".
#' @keywords internal
#' @noRd
alignSpectra_landmark_cpp <- function(spectra_cal,
                                      peaks_all,
                                      ref_peaks,
                                      tolerance_align = 0.002,
                                      ppm = FALSE,
                                      n_workers = max(1, parallel::detectCores() - 1),
                                      verbose = TRUE) {
  n_workers <- .safe_n_workers(n_workers, fallback = 1L)
  stopifnot(length(spectra_cal) == length(peaks_all))
  if (!exists("cpp_match_anchors", mode = "function") ||
      !exists("cpp_warp_mz_piecewise", mode = "function")) {
    stop("C++ kernels cpp_match_anchors/cpp_warp_mz_piecewise not found. ",
         "Ensure src/landmark_align.cpp is compiled and exported with Rcpp.")
  }
  if (verbose) message(sprintf("[alignSpectra_landmark_cpp] Aligning %d spectra", length(spectra_cal)))

  ref_mz <- MALDIquant::mass(ref_peaks)

  worker_fun <- function(idx) {
    s  <- spectra_cal[[idx]]
    pk <- peaks_all[[idx]]

    src <- MALDIquant::mass(pk)
    if (length(src)) src <- sort(src)

    mt  <- cpp_match_anchors(src, ref_mz, tol = tolerance_align, ppm = ppm)
    mz  <- MALDIquant::mass(s)
    iz  <- MALDIquant::intensity(s)
    mz2 <- cpp_warp_mz_piecewise(mz, mt$src, mt$dst)

    MALDIquant::createMassSpectrum(
      mass = mz2,
      intensity = iz,
      metaData = MALDIquant::metaData(s)
    )
  }

  ids <- seq_along(spectra_cal)
  if (n_workers <= 1L) return(lapply(ids, worker_fun))

  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE))
  parallel::clusterEvalQ(cl, { library(MALDIquant) })
  parallel::clusterExport(
    cl,
    varlist = c("spectra_cal", "peaks_all", "ref_mz", "tolerance_align", "ppm"),
    envir = environment()
  )
  parallel::parLapply(cl, ids, worker_fun)
}
