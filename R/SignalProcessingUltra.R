#' Optimized signal processing for MALDI-TOF spectra (parallel + optional C++ alignment)
#'
#' This function performs post-acquisition processing for lists of MALDIquant MassSpectrum objects.
#' It parallelizes per-spectrum steps and can align with either MALDIquant's warping (cubic/lowess)
#' or a fast landmark-based C++ algorithm ("landmark_cpp").
#'
#' @param x list of MALDIquant MassSpectrum objects.
#' @param transformIntensity_method character, intensity transform (default "log").
#' @param smoothing_method character, smoothing method ("Wavelet" UDWT).
#' @param removeBaseline_method character, baseline method ("TopHat" default; "SNIP","ConvexHull" supported).
#' @param removeBaseline_iterations integer, SNIP iterations if removeBaseline_method = "SNIP".
#' @param calibrateIntensity_method character, intensity calibration ("PQN" default, or "TIC","median").
#' @param alignSpectra_NoiseMethod character, noise estimator for peak finding pre-alignment ("MAD").
#' @param alignSpectra_method character, alignment engine: "cubic" (default), "lowess", or "landmark_cpp".
#' @param alignSpectra_halfWs integer, half window size for peak detection.
#' @param alignSpectra_SN numeric, SNR for peak detection.
#' @param tolerance_align numeric, tolerance for matching anchors to the reference during alignment.
#'                       Use consistent units across your pipeline (Da by default here).
#' @param ppm_align logical, set TRUE if tolerance_align is in ppm (then interpreted as ppm).
#' @param referenceSpectra optional MALDIquant MassPeaks object to use as alignment reference.
#' @param minFrequency numeric, minimum peak frequency to build reference if not provided (default 0.7).
#' @param binPeaks_method character, "strict" (default) or "relaxed" for reference peak binning.
#' @param keepReferenceSpectra logical, if TRUE and no reference provided, returns list(spectra=..., RefS=...).
#' @param n_workers integer, number of parallel workers (default: all cores minus one).
#' @param ref_sample_n integer or NULL, if set, build the reference from a random subset of this many spectra.
#' @param verbose logical, print progress.
#' @param ... passed to MALDIrppa::wavSmoothing (e.g., n.levels).
#'
#' @return A list of MassSpectrum objects, or a list with $spectra and $RefS if keepReferenceSpectra = TRUE.
#' @export
SignalProcessingUltra <- function(
    x,
    transformIntensity_method = "log",
    smoothing_method = "Wavelet",
    removeBaseline_method = "TopHat",
    removeBaseline_iterations = 25,
    calibrateIntensity_method = "PQN",
    alignSpectra_NoiseMethod = "MAD",
    alignSpectra_method = c("cubic", "lowess", "landmark_cpp"),
    alignSpectra_halfWs = 11,
    alignSpectra_SN = 3,
    tolerance_align = 0.002,   # keep unit consistent across your pipeline (Da by default)
    ppm_align = FALSE,         # set TRUE if tolerance_align is in ppm
    referenceSpectra = NULL,
    minFrequency = 0.7,
    binPeaks_method = "strict",
    keepReferenceSpectra = FALSE,
    n_workers = NULL,          # explicit cores for Unix; ignored on Windows (stays serial)
    ref_sample_n = NULL,       # e.g., 2000 to speed up reference building
    verbose = TRUE,
    ...
) {
  stopifnot(is.list(x))
  alignSpectra_method <- match.arg(alignSpectra_method)
  n_workers <- .safe_n_workers(n_workers, fallback = 1L)

  is_windows <- .Platform$OS.type == "windows"
  os_label <- if (is_windows) "windows" else "unix"
  # For speed and stability on Windows, keep serial/vectorized path
  can_parallel_unix <- (!is_windows) && (n_workers > 1L)

  if (verbose) {
    message(sprintf(
      "[SignalProcessingUltra] N=%d; workers=%d; align=%s; OS=%s",
      length(x), as.integer(n_workers), alignSpectra_method, os_label
    ))
  }

  # Convenience: mclapply only on Unix; else fall back to lapply
  u_mclapply <- function(obj, fun, cores) {
    if (can_parallel_unix) {
      parallel::mclapply(obj, fun, mc.cores = n_workers)
    } else {
      lapply(obj, fun)
    }
  }

  # 1) Transform (batch)
  if (verbose) message("1) Transform: ", transformIntensity_method)
  spectra <- MALDIquant::transformIntensity(x, method = transformIntensity_method)

  # 2) Smooth (batch; MALDIrppa accepts a list of MassSpectrum)
  if (verbose) message("2) Smooth: ", smoothing_method)
  spectra <- MALDIrppa::wavSmoothing(spectra, method = smoothing_method, ...)

  # 3) Baseline (batch)
  if (verbose) message("3) Baseline: ", removeBaseline_method)
  spectra <- switch(removeBaseline_method,
                    TopHat     = MALDIquant::removeBaseline(spectra, method = "TopHat"),
                    SNIP       = MALDIquant::removeBaseline(spectra, method = "SNIP", iterations = removeBaseline_iterations),
                    ConvexHull = MALDIquant::removeBaseline(spectra, method = "ConvexHull"),
                    stop("Unsupported removeBaseline_method: ", removeBaseline_method)
  )

  # 4) Calibrate (PQN on list; TIC/median also on list)
  if (identical(calibrateIntensity_method, "PQN")) {
    if (verbose) message("4) Calibrate: PQN (batch)")
    spectra_cal <- if (length(spectra) < 2L) {
      warning("PQN not supported for a single spectrum; falling back to TIC.")
      MALDIquant::calibrateIntensity(spectra, method = "TIC")
    } else {
      MALDIquant::calibrateIntensity(spectra, method = "PQN")
    }
  } else {
    if (verbose) message("4) Calibrate: ", calibrateIntensity_method)
    spectra_cal <- MALDIquant::calibrateIntensity(spectra, method = calibrateIntensity_method)
  }

  # 5) Build reference peaks (from subset for speed if requested)
  if (is.null(referenceSpectra)) {
    if (verbose) message("5) Building reference peaks (minFrequency=", minFrequency, ", method=", binPeaks_method, ")")
    idx_ref <- seq_along(spectra_cal)
    if (!is.null(ref_sample_n) && ref_sample_n < length(spectra_cal)) {
      set.seed(1L)
      idx_ref <- sort(sample.int(length(spectra_cal), ref_sample_n))
      if (verbose) message("   Using ", length(idx_ref), " spectra to derive reference")
    }
    peaks_ref <- MALDIquant::detectPeaks(
      spectra_cal[idx_ref],
      halfWindowSize = alignSpectra_halfWs,
      method = alignSpectra_NoiseMethod,
      SNR = alignSpectra_SN
    )
    ref_peaks <- MALDIquant::referencePeaks(peaks_ref, minFrequency = minFrequency, method = binPeaks_method)
  } else {
    if (verbose) message("5) Using provided reference peaks")
    ref_peaks <- referenceSpectra
  }

  # 6) Align
  if (verbose) message("6) Align: ", alignSpectra_method)
  if (alignSpectra_method %in% c("cubic", "lowess")) {
    # MALDIquant alignment (batch, fastest/most stable on Windows)
    aligned <- MALDIquant::alignSpectra(
      spectra_cal,
      halfWindowSize  = alignSpectra_halfWs,
      SNR             = alignSpectra_SN,
      tolerance       = tolerance_align,
      warpingMethod   = alignSpectra_method,
      allowNoMatches  = TRUE,
      emptyNoMatches  = FALSE,
      reference       = ref_peaks
    )
  } else {
    # Fast landmark alignment via Rcpp; compute anchors in batch
    peaks_all <- MALDIquant::detectPeaks(
      spectra_cal,
      halfWindowSize = alignSpectra_halfWs,
      method = alignSpectra_NoiseMethod,
      SNR = alignSpectra_SN
    )
    aligned <- alignSpectra_landmark_cpp(
      spectra_cal     = spectra_cal,
      peaks_all       = peaks_all,
      ref_peaks       = ref_peaks,
      tolerance_align = tolerance_align,
      ppm             = ppm_align,
      # Keep serial on Windows (to avoid PSOCK); allow multiple cores on Unix
      n_workers       = if (is_windows) 1L else n_workers,
      verbose         = verbose
    )
  }

  if (keepReferenceSpectra && is.null(referenceSpectra)) {
    return(list(spectra = aligned, RefS = ref_peaks))
  }
  aligned
}
