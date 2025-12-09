#' Signal processing for MALDI-TOF spectra (wrapper to SignalProcessingUltra)
#'
#' Backward-compatible wrapper that delegates to SignalProcessingUltra.
#' Keeps the original argument names/signature so existing code continues to work.
#'
#' @inheritParams SignalProcessingUltra
#' @param ... additional arguments passed to SignalProcessingUltra (e.g., n_workers, ref_sample_n).
#' @return A list of processed MassSpectrum objects, or list(spectra, RefS) if keepReferenceSpectra = TRUE.
#' @export
SignalProcessing <- function(
    x,
    transformIntensity_method = "sqrt",
    smoothing_method = "Wavelet",
    removeBaseline_method = "SNIP",
    removeBaseline_iterations = 25,
    calibrateIntensity_method = "TIC",
    alignSpectra_NoiseMethod = "MAD",
    alignSpectra_method = "lowess",
    alignSpectra_halfWs = 11,
    alignSpectra_SN = 3,
    tolerance_align = 0.002,
    referenceSpectra = NULL,
    minFrequency = 0.5,
    binPeaks_method = "strict",
    keepReferenceSpectra = FALSE,
    ...
) {
  # Map to Ultra; accept extra speed knobs via ...
  SignalProcessingUltra(
    x                           = x,
    transformIntensity_method   = transformIntensity_method,
    smoothing_method            = smoothing_method,
    removeBaseline_method       = removeBaseline_method,
    removeBaseline_iterations   = removeBaseline_iterations,
    calibrateIntensity_method   = calibrateIntensity_method,
    alignSpectra_NoiseMethod    = alignSpectra_NoiseMethod,
    alignSpectra_method         = alignSpectra_method,
    alignSpectra_halfWs         = alignSpectra_halfWs,
    alignSpectra_SN             = alignSpectra_SN,
    tolerance_align             = tolerance_align,
    referenceSpectra            = referenceSpectra,
    minFrequency                = minFrequency,
    binPeaks_method             = binPeaks_method,
    keepReferenceSpectra        = keepReferenceSpectra,
    ...
  )
}