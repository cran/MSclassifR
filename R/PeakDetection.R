#' Detection of peaks in MassSpectrum objects
#'
#' Detects peaks on a list of MALDIquant MassSpectrum objects, with an optional
#' preliminary averaging step and an optional discrete-bin alignment of detected
#' peaks. Per-spectrum peak detection can be parallelized on Unix-alike systems
#' (Linux/macOS) for very large inputs; on Windows a serial/vectorized path is
#' used. The result is a list of MALDIquant MassPeaks ready for downstream
#' matrix building (e.g., with MALDIquant::intensityMatrix or build_X_from_peaks_fast).
#'
#' @param x List of MALDIquant::MassSpectrum objects (one per sample). These are
#'   typically obtained after preprocessing (baseline, smoothing, normalization).
#' @param averageMassSpec Logical; if TRUE, average spectra using
#'   MALDIquant::averageMassSpectra before peak detection. If `labels` is provided
#'   and its length equals length(x), a groupwise averaging is performed; otherwise
#'   all spectra are averaged. Default TRUE.
#' @param labels Optional factor/character vector for groupwise averaging (same
#'   semantics as MALDIquant::averageMassSpectra). Ignored if `averageMassSpec = FALSE`.
#' @param averageMassSpectraMethod Character, "median" (default) or "mean". Passed
#'   to MALDIquant::averageMassSpectra when `averageMassSpec = TRUE`.
#' @param SNRdetection Numeric; signal-to-noise ratio threshold for peak detection
#'   (MALDIquant::detectPeaks argument SNR). Default 3.
#' @param binPeaks Logical; if TRUE, align detected peaks into discrete bins using
#'   MALDIquant::binPeaks (method/tolerance set by `AlignMethod`/`Tolerance`). Default TRUE.
#' @param PeakDetectionMethod Character; MALDIquant::detectPeaks method, e.g.,
#'   "MAD" (default) or "SuperSmoother".
#' @param halfWindowSizeDetection Integer; half window size for local maxima
#'   (MALDIquant::detectPeaks argument halfWindowSize). Default 11.
#' @param AlignMethod Character; MALDIquant::binPeaks method, "strict" (default)
#'   or "relaxed".
#' @param Tolerance Numeric; MALDIquant::binPeaks tolerance (units consistent with
#'   your m/z axis). Default 0.002.
#' @param n_workers Integer or NULL; requested number of parallel workers for the
#'   Unix mclapply path. On Windows, this is ignored (serial path). The effective
#'   number is sanitized by an internal helper (.safe_n_workers) to avoid
#'   oversubscription and R CMD check issues. Default NULL (auto).
#' @param verbose Logical; if TRUE, print progress messages. Default TRUE.
#' @param min_parallel_n Integer; minimum number of spectra at which the function
#'   will attempt Unix parallelization via mclapply. Defaults to 2000. Increase to
#'   be more conservative, decrease to parallelize more aggressively. Set to Inf
#'   to effectively disable Unix parallelization regardless of n_workers.
#' @param chunk_size Integer; number of spectra per chunk/task submitted to
#'   mclapply on Unix. Larger chunks reduce scheduling overhead but use more
#'   memory per task. Default 1000.
#' @param ... Reserved for future extensions or pass-through to MALDIquant/MALDIrppa.
#'
#' @return A list of MALDIquant::MassPeaks objects (one per input or averaged
#'   spectrum). If `binPeaks = TRUE`, all MassPeaks are aligned to the same
#'   discrete m/z bins (shared centers), facilitating fast matrix construction.
#'
#' @details
#' - Averaging: if `averageMassSpec = TRUE` and `labels` is provided with
#'   length(labels) == length(x), MALDIquant::averageMassSpectra performs a
#'   groupwise averaging by labels. Otherwise, all spectra are averaged. If
#'   `averageMassSpec = FALSE`, the input list is used as-is.
#' - Peak detection: peak finding uses MALDIquant::detectPeaks with the given
#'   SNR/method/half-window. This is applied per spectrum (serial or parallel).
#' - Discrete-bin alignment: when `binPeaks = TRUE`, MALDIquant::binPeaks aligns
#'   detected peaks to a shared discrete grid (method `AlignMethod`, tolerance
#'   `Tolerance`), enabling consistent feature columns across spectra.
#' - Parallelization: on Windows, a single serial/vectorized call to
#'   MALDIquant::detectPeaks is used (fast enough for small/medium inputs). On
#'   Unix-alike systems, when `length(x) >= min_parallel_n` and `n_workers > 1`,
#'   the list is split into chunks of size `chunk_size` and processed with
#'   parallel::mclapply using the requested number of workers.
#' - Meta-data: if `labels` is provided, label information is appended best-effort
#'   to MassPeaks metaData (file/fullName), preserving existing fields where possible.
#'
#' @examples
#' if (requireNamespace("MALDIquant", quietly = TRUE)) {
#'   # Two toy spectra with peaks near 1000, 1500, 2000 Da
#'   mass <- seq(900, 2100, by = 1)
#'   make_spectrum <- function(shift) {
#'     inten <- dnorm(mass, 1000 + shift, 2) * 50 +
#'              dnorm(mass, 1500 - shift, 2) * 80 +
#'              dnorm(mass, 2000 + shift, 2) * 40 +
#'              rnorm(length(mass), 0, 0.2)
#'     MALDIquant::createMassSpectrum(mass = mass, intensity = inten)
#'   }
#'   spectra <- list(make_spectrum(0.3), make_spectrum(-0.3))
#'
#'   # Detect peaks without averaging; align in strict bins
#'   peaks <- PeakDetection(
#'     x = spectra,
#'     averageMassSpec = FALSE,
#'     SNRdetection = 3,
#'     PeakDetectionMethod = "MAD",
#'     binPeaks = TRUE,
#'     AlignMethod = "strict",
#'     Tolerance = 0.5,
#'     verbose = TRUE
#'   )
#'
#'   # Build an intensity matrix (rows = spectra, cols = aligned m/z bins)
#'   X <- MALDIquant::intensityMatrix(peaks)
#'   dim(X)
#' }
#'
#' @seealso MALDIquant::averageMassSpectra, MALDIquant::detectPeaks,
#'   MALDIquant::binPeaks; MALDIquant::intensityMatrix; build_X_from_peaks_fast
#'   for a fast matrix builder from MassPeaks.
#' @export
PeakDetection <- function(
    x,
    averageMassSpec = TRUE,
    labels = NULL,
    averageMassSpectraMethod = "median",
    SNRdetection = 3,
    binPeaks = TRUE,
    PeakDetectionMethod = "MAD",
    halfWindowSizeDetection = 11,
    AlignMethod = "strict",
    Tolerance = 0.002,
    n_workers = NULL,            # used only on Unix; ignored on Windows
    verbose = TRUE,
    min_parallel_n = 2000L,      # only consider Unix parallel when very large
    chunk_size = 1000L,
    ...
) {
  stopifnot(is.list(x))
  n_workers <- .safe_n_workers(n_workers, fallback = 1L)
  is_windows <- .Platform$OS.type == "windows"
  can_parallel_unix <- (!is_windows) && (n_workers > 1L)

  # Averaging
  if (isTRUE(averageMassSpec)) {
    if (verbose) message("Average of the MassSpectrum objects")
    if (!is.null(labels) && length(labels) == length(x)) {
      message("The length of 'labels' equals the number of spectra. Consider averageMassSpec = FALSE if you do not want averaging.")
    }
    avgSpectra <- MALDIquant::averageMassSpectra(x, labels = labels, method = averageMassSpectraMethod)
  } else {
    avgSpectra <- x
  }

  n_spec <- length(avgSpectra)
  # Decide strategy: Windows -> serial; Unix -> serial unless huge
  use_unix_parallel <- can_parallel_unix && (n_spec >= min_parallel_n)

  if (verbose) {
    message(sprintf("Detecting peaks (method=%s, SNR=%s) [%s]",
                    PeakDetectionMethod, SNRdetection,
                    if (use_unix_parallel) sprintf("parallel mclapply x%d", n_workers) else "serial"))
  }

  if (use_unix_parallel) {
    # Chunk indices and process with mclapply (no PSOCK serialization)
    idx <- split(seq_len(n_spec), ceiling(seq_len(n_spec) / as.integer(chunk_size)))
    chunks <- parallel::mclapply(idx, function(ii) {
      MALDIquant::detectPeaks(
        avgSpectra[ii],
        SNR = SNRdetection,
        method = PeakDetectionMethod,
        halfWindowSize = halfWindowSizeDetection
      )
    }, mc.cores = n_workers)
    peaks <- do.call(c, chunks)
  } else {
    # Serial, batch call (fastest on Windows and small/medium jobs)
    peaks <- MALDIquant::detectPeaks(
      avgSpectra,
      SNR = SNRdetection,
      method = PeakDetectionMethod,
      halfWindowSize = halfWindowSizeDetection
    )
  }

  if (isTRUE(binPeaks)) {
    if (verbose) message("Aligning peaks in discrete bins")
    peaks <- MALDIquant::binPeaks(peaks, method = AlignMethod, tolerance = Tolerance)
  }

  # Add labels to metaData if provided (best-effort)
  if (!is.null(labels)) {
    nlab <- length(labels); npk <- length(peaks)
    for (i in seq_len(npk)) {
      lb <- if (nlab >= i) labels[[i]] else NA_character_
      if (is.null(MALDIquant::metaData(peaks[[i]])$fullName)) {
        peaks[[i]]@metaData$file <- paste0(lb, ".", MALDIquant::metaData(peaks[[i]])$file %||% "")
      } else {
        peaks[[i]]@metaData$fullName <- paste0(lb, ".", MALDIquant::metaData(peaks[[i]])$fullName)
      }
    }
  }

  peaks
}

# small infix helper to handle NULLs in paste
`%||%` <- function(a, b) if (is.null(a)) b else a
