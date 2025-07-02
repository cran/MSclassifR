
SignalProcessing<-function(x,
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
                           ...) {
  ## Signal processing ##

  # Transform intensity
  print(paste0("1- Transform intensities with ", transformIntensity_method,
               " method"))
  spectra <- MALDIquant::transformIntensity(x, method = transformIntensity_method)

  # Smoothing
  print(paste0("2- Smoothing with ", smoothing_method, " method"))
  spectra <- MALDIrppa::wavSmoothing(spectra, method = smoothing_method,...)

  # Remove Baseline
  print(paste0("3- Remove baseline with ", removeBaseline_method,
               " method"))
  if(removeBaseline_method == "SNIP"){
    spectra <- MALDIquant::removeBaseline(spectra, method = removeBaseline_method,
                                          iterations = removeBaseline_iterations);}
  if(removeBaseline_method == "TopHat"){
    spectra <- MALDIquant::removeBaseline(spectra, method = removeBaseline_method);}

  # Calibration intensities
  print(paste0("4- Calibrate intensity with ", calibrateIntensity_method,
               " method"))
  spectra_before <- MALDIquant::calibrateIntensity(spectra,
                                                   method = calibrateIntensity_method)

  #Reference mass spectrum for alignment
  peaks_before <- MALDIquant::detectPeaks(spectra_before,
                                          halfWindowSize = alignSpectra_halfWs,
                                          method = alignSpectra_NoiseMethod,
                                          SNR = alignSpectra_SN)
  if (is.null(referenceSpectra)) {
    print(paste0("5- Creating a reference spectrum to align spectra"))
    referenceSpectraAlign <- MALDIquant::referencePeaks(peaks_before,
                                                        minFrequency = minFrequency, method = binPeaks_method)
  }
  else {
    print(paste0("5- Reference spectrum to align the spectra has been provided"))
    referenceSpectraAlign <- referenceSpectra
  }

  # Alignment
  print(paste0("6- Align spectra with ", alignSpectra_method,
               " method"))
  spectraPro <- try(MALDIquant::alignSpectra(spectra_before,
                                             halfWindowSize = alignSpectra_halfWs,
                                             SNR = alignSpectra_SN,
                                             tolerance = tolerance_align,
                                             warpingMethod = alignSpectra_method,
                                             allowNoMatches = TRUE,
                                             emptyNoMatches = FALSE,
                                             reference = referenceSpectraAlign))

  if (is.character(spectraPro)) {
    warning("Processed spectra are available but were not aligned due to an error (minFrequency is maybe too high)!")
    spectra <- spectra_before
  }

  if (is.null(referenceSpectra)){
    if (keepReferenceSpectra == TRUE){
       spectra = list(spectra = spectraPro, RefS = referenceSpectraAlign);
    }else {spectra <- spectraPro;}
  }else {spectra <- spectraPro;}

  return(spectra)
}
