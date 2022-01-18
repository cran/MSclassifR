\name{PeakDetection}

\alias{PeakDetection}

\title{Detection of peaks in \code{MassSpectrum} objects.}

\description{
This function performs a data analysis pipeline to pre-process mass spectra. It provides average intensities and detects peaks using functions of R packages \code{MALDIquant} and \code{MALDIrppa}.
}

\usage{
PeakDetection(x,
              labels,
              averageMassSpectraMethod = "median",
              SNRdetection = 3,
              halfWindowSizeDetection = 11,
              AlignFrequency = 0.20,
              AlignMethod = "strict",
              Tolerance = 0.002,
              ...)

}

\arguments{
  \item{x}{a \code{list} of \code{MassSpectrum} objects (see \code{MALDIquant} R package).}

  \item{labels}{a \code{list} of \code{factor} objects to do groupwise averaging.}

  \item{averageMassSpectraMethod}{a \code{character} indicating the method used to average mass spectra according to labels. It is fixed to \code{"median"} by default. See \code{averageMassSpectra} of \code{MALDIquant} R package.}

  \item{SNRdetection}{a \code{numeric} value indicating the signal-to-noise ratio used to detect peaks (by default = {\code{3}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{halfWindowSizeDetection}{a \code{numeric} value half window size to detect peaks (by default = {\code{11}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{AlignFrequency}{a \code{numeric} value used to align and bin mass spectra using minimum relative frequency. See \code{alignPeaks} of the \code{MALDIrppa} R package for more details.}

  \item{AlignMethod}{a \code{character} indicating the method used to equalize masses for similar peaks. The {\code{strict}} method is used by default. See \code{binPeaks} of the \code{MALDIquant} R package for more details.}

  \item{Tolerance}{a \code{numeric} value corresponding to the maximal deviation in peak masses to be considered as identical (by default = {\code{0.002}}). See \code{determineWarpingFunctions} of the \code{MALDIquant} R package for details.}

  \item{...}{other arguments from \code{MALDIquant} and \code{MALDIrppa} packages.}
}

\value{Returns a list of \code{MassPeaks} objects (see \code{MALDIquant} R package) for each mass spectrum in \code{x}.}

\references{
Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. Bioinformatics. 2012 Sep 1;28(17):2270-1. \doi{10.1093/bioinformatics/bts447}. Epub 2012 Jul 12. PMID: 22796955.

Javier Palarea-Albaladejo, Kevin Mclean, Frank Wright, David G E Smith, MALDIrppa: quality control and robust analysis for mass spectrometry data, Bioinformatics, Volume 34, Issue 3, 01 February 2018, Pages 522 - 523, \doi{https://doi.org/10.1093/bioinformatics/btx628}
}

\seealso{
\bold{Vignette MSclassifR :}
\url{https://agodmer.github.io/StageDubii2021/Package_MSclassifR-in-progress/MSclassifR/vignettes/Vignette-MSclassifR.html}
}

\examples{
library("MALDIquant")
library("MALDIquantForeign")
library("MSclassifR")


## Load mass spectra and metadata
data("CitrobacterRKIspectra", "CitrobacterRKImetadata", package = "MSclassifR")

## Pre-processing of mass spectra
spectra <- SignalProcessing(CitrobacterRKIspectra)

## Detection of peaks in pre-processed mass spectra
peaks <- PeakDetection(x = spectra,
 labels = CitrobacterRKImetadata$Strain_name_spot,
 averageMassSpectraMethod = "median",
 SNRdetection = 3,
 halfWindowSizeDetection = 11,
 AlignFrequency = 0.20,
 AlignMethod = "strict",
 Tolerance = 0.002)

# Plot peaks on a pre-processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],col_spec="blue",col_peak="black")
}

