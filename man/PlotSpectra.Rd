\name{PlotSpectra}

\alias{PlotSpectra}

\title{Plot mass spectra with detected peaks}

\description{
This function performs a plot of a \code{AbstractMassObject} object (see the \code{MALDIquant} R package). It can be used to highlight peaks in a mass spectrum.}


\usage{

PlotSpectra(SpectralData, absx="ALL", Peaks=NULL, Peaks2=NULL, col_spec=1, 
            col_peak=2, shape_peak=3, col_peak2=2, shape_peak2=2)

}

\arguments{
  \item{SpectralData}{\code{MassSpectrum} object of S4 class (see \code{MALDIquant} R package).}
  \item{absx}{\code{vector} indicating lower and upper bounds for the mass-over-charge values to plot.}
  \item{Peaks}{\code{MassPeaks} object (see \code{MALDIquant} R package). If NULL, peaks are not highlighted.}
  \item{Peaks2}{numeric \code{vector} of mass-over-charge values to plot on the mass spectrum.}
  \item{col_spec}{color of the mass spectrum.}
  \item{col_peak}{color of the peak points corresponding to \code{Peaks}.}
  \item{shape_peak}{shape of the peak points corresponding to \code{Peaks}.}
  \item{col_peak2}{color of the peak points corresponding to \code{Peaks2}.}
  \item{shape_peak2}{Shape of the peak points corresponding to \code{Peaks2}.}
}

\value{A \code{ggplot} object (see \code{ggplot2} R package). Mass-over-charge values are in x-axis and intensities in y-axis.}


\examples{

library("MSclassifR")

# Load mass spectra
data("CitrobacterRKIspectra", package = "MSclassifR")
# Plot raw mass spectrum
PlotSpectra(SpectralData = CitrobacterRKIspectra[[1]])
# standard pre-processing of mass spectra
spectra <- SignalProcessing(CitrobacterRKIspectra)
# Plot pre-processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]])
# detection of peaks in pre-processed mass spectra
peaks <- PeakDetection(x = spectra, averageMassSpec=FALSE)
# Plot peaks on pre-processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],col_spec="blue",col_peak="black")

}



