## Peak Detection ##

PeakDetection <- function(x,
                          averageMassSpec=TRUE,
                          labels=NULL,
                          averageMassSpectraMethod="median",
                          SNRdetection=3,
                          binPeaks = TRUE,
                          PeakDetectionMethod="MAD",
                          halfWindowSizeDetection=11,
                          AlignMethod="strict",
                          Tolerance=0.002,
                          ...)

{

  ## Peaks detection =================================
  #  Average spectra
  if(averageMassSpec==TRUE){
  print("Average of the MassSpectrum objects")
  if(length(labels)==length(x)){message("The length of the labels argument is equal to the length of the spectra, it is advisable to replace averageMassSpec=TRUE with averageMassSpec=FALSE")}
  avgSpectra <- MALDIquant::averageMassSpectra(x,labels = labels ,method = averageMassSpectraMethod)
  }else{avgSpectra=x}
  # Peak detection
  peaks <- MALDIquant::detectPeaks(avgSpectra, SNR= SNRdetection, method = PeakDetectionMethod, halfWindowSize = halfWindowSizeDetection)

  ## Peaks alignment and bin  =================================
  # Alignment peaks
  #peaks <- MALDIrppa::alignPeaks(peaks, minFreq = AlignFrequency, tolerance = Tolerance)
  # Collect peaks
  if(binPeaks==TRUE){
  print("Aligning peaks in discrete bins")
  peaks <- MALDIquant::binPeaks(peaks, method = AlignMethod, tolerance = Tolerance)
  }else{peaks=peaks}
  
  #renaming fullName using labels
  for (i in 1:length(peaks)){
  if(is.null(peaks[[i]]@metaData$fullName)){peaks[[i]]@metaData$file=paste(c(labels[i],".",peaks[[i]]@metaData$file),collapse="")} 
  else{peaks[[i]]@metaData$fullName=paste(c(labels[i],".",peaks[[i]]@metaData$fullName),collapse="")}
  }

  return(peaks)
}
