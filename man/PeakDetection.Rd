\name{PeakDetection}

\alias{PeakDetection}

\title{Detection of peaks in \code{MassSpectrum} objects.}

\description{
This function performs a data analysis pipeline to pre-process mass spectra. It provides average intensities and detects peaks using functions of R packages \code{MALDIquant} and \code{MALDIrppa}.
}

\usage{
PeakDetection(x,
              averageMassSpec = TRUE,
              labels = NULL,
              averageMassSpectraMethod = "median",
              SNRdetection = 3,
              binPeaks = TRUE,
              PeakDetectionMethod = "MAD",
              halfWindowSizeDetection = 11,
              AlignMethod = "strict",
              Tolerance = 0.002,
              ...)

}

\arguments{
  \item{x}{a \code{list} of \code{MassSpectrum} objects (see \code{MALDIquant} R package).}

  \item{averageMassSpec}{a \code{logical} value indicating whether it is necessary to group the spectra according to the arguments \code{labels} and  \code{averageMassSpectraMethod}. It is fixed to \code{"TRUE"} by default.}

  \item{labels}{a \code{list} of \code{factor} objects to do groupwise averaging.}

  \item{averageMassSpectraMethod}{a \code{character} indicating the method used to average mass spectra according to labels. It is fixed to \code{"median"} by default.This function can be replaced by another mathematical function such as \code{"mean"}. See \code{averageMassSpectra} of \code{MALDIquant} R package.}

  \item{PeakDetectionMethod}{a \code{character} indicating the noise estimation method. It uses \code{"MAD"} method for \code{list} of \code{MassSpectrum} objects. This noise estimation method estimation method can be remplaced \code{"SuperSmoother"}. See \code{estimateNoise-methods} of the \code{MALDIquant} R package for details.}
  
  \item{SNRdetection}{a \code{numeric} value indicating the signal-to-noise ratio used to detect peaks (by default = {\code{3}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}
  
  \item{binPeaks}{a \code{logical} value indicating the peaks are aligned in discrete bins. It is fixed to \code{"TRUE"} by default. See \code{binPeaks} of the \code{MALDIquant} R package for details.}

  \item{halfWindowSizeDetection}{a \code{numeric} value half window size to detect peaks (by default = {\code{11}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{AlignMethod}{a \code{character} indicating the method used to equalize masses for similar peaks. The \code{"strict"} method is used by default corresponding to a unique peak per bin from the same sample. This method can be remplaced by \code{"relaxed"} corresponding to multiple peaks per bin from the same sample. See \code{binPeaks} of the \code{MALDIquant} R package for more details.}

  \item{Tolerance}{a \code{numeric} value corresponding to the maximal deviation in peak masses to be considered as identical in ppm (by default = {\code{0.002}}). See \code{determineWarpingFunctions} of the \code{MALDIquant} R package for details.}

  \item{...}{other arguments from \code{MALDIquant} and \code{MALDIrppa} packages.}
}

\value{Returns a list of \code{MassPeaks} objects (see \code{MALDIquant} R package) for each mass spectrum in \code{x}.}

\references{
Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. Bioinformatics. 2012 Sep 1;28(17):2270-1. \doi{10.1093/bioinformatics/bts447}. Epub 2012 Jul 12. PMID: 22796955.

Javier Palarea-Albaladejo, Kevin Mclean, Frank Wright, David G E Smith, MALDIrppa: quality control and robust analysis for mass spectrometry data, Bioinformatics, Volume 34, Issue 3, 01 February 2018, Pages 522 - 523, \doi{https://doi.org/10.1093/bioinformatics/btx628}
}

\details{
The \code{PeakDetection} function provides an analysis pipeline for \code{MassSpectrum} objects including peaks detection and binning.

All the methods used for \code{PeakDetection} functions are selected from \code{MALDIquant} and \code{MALDIrppa} packages.
}


\examples{
library("MALDIquant")
library("MSclassifR")


## Load mass spectra and metadata
data("CitrobacterRKIspectra", "CitrobacterRKImetadata", package = "MSclassifR")

## Pre-processing of mass spectra
spectra <- SignalProcessing(CitrobacterRKIspectra)

## Detection of peaks in pre-processed mass spectra
peaks <- PeakDetection(x = spectra,
                       averageMassSpec = FALSE,
                       labels = CitrobacterRKImetadata$Strain_name_spot,
                       averageMassSpectraMethod = "median",
                       SNRdetection = 3,
                       binPeaks = TRUE,
                       halfWindowSizeDetection = 11,
                       AlignFrequency = 0.20,
                       AlignMethod = "strict",
                       Tolerance = 0.002)

# Plot peaks on a pre-processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],col_spec="blue",col_peak="black")
}

\keyword{Preprocessing}
