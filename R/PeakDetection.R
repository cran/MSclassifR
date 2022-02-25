## Peak Detection ##

PeakDetection <- function(x,
                          labels,
                          averageMassSpectraMethod = "median",
                          SNRdetection = 3,
                          halfWindowSizeDetection = 11,
                          AlignFrequency = 0.20,
                          AlignMethod = "strict",
                          Tolerance = 0.002,
                          ...)

{

  ## Peaks detection =================================
  #  Average spectra
  avgSpectra <- MALDIquant::averageMassSpectra(x,labels = labels ,method = averageMassSpectraMethod )
  # Peak detection
  peaks <- MALDIquant::detectPeaks(avgSpectra, SNR= SNRdetection, halfWindowSize = halfWindowSizeDetection)

  ## Peaks alignment and bin  =================================
  # Alignment peaks
  peaks <- MALDIrppa::alignPeaks(peaks, minFreq = AlignFrequency, tolerance = Tolerance)
  # Collect peaks
  peaks <- MALDIquant::binPeaks(peaks, method=c("strict"), tolerance = Tolerance)

  #renaming fullName using labels
  for (i in 1:length(peaks)){
  peaks[[i]]@metaData$fullName=paste(c(labels[i],".",peaks[[i]]@metaData$fullName),collapse="")
  }

  return(peaks)
}
