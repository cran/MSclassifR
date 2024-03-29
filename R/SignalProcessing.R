########################
## Spectral treatment ##
########################

SignalProcessing <- function (x,
                              transformIntensity_method = "sqrt",
                              smoothing_method = "Wavelet",
                              removeBaseline_method = "SNIP",
                              removeBaseline_iterations = 25,
                              calibrateIntensity_method = "TIC",
                              alignSpectra_halfWs = 11,
                              alignSpectra_SN = 3,
                              tolerance_align = 0.002,
                              ...){
  ## Signal treatment ##

  # Transform intensity
  spectra <- MALDIquant::transformIntensity(x,
                                            method = transformIntensity_method)
  # Smoothing
  spectra <- MALDIrppa::wavSmoothing(spectra,
                                     method = smoothing_method)
  # Remove Baseline
  spectra <- MALDIquant::removeBaseline(spectra,
                                        method = removeBaseline_method,
                                        iterations = removeBaseline_iterations)
  # Calibration intensities
  spectra <- MALDIquant::calibrateIntensity(spectra,
                                            method = calibrateIntensity_method)
  #Alignment
  spectra <- MALDIquant::alignSpectra(spectra,
                                      halfWindowSize = alignSpectra_halfWs,
                                      SNR = alignSpectra_SN, tolerance = tolerance_align,
                                      allowNoMatches = TRUE,
                                      emptyNoMatches = FALSE)

  return(spectra) # List of treated spectra
}
