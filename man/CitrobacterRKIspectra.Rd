\docType{data}

\name{CitrobacterRKIspectra}

\alias{CitrobacterRKIspectra}

\title{Mass spectra corresponding to the bacterial species \emph{Citrobacter} sp. from  The Robert Koch-Institute (RKI) database of microbial MALDI-TOF mass spectra}

\format{
A list that contains 14 objects of class S4 corresponding each to a each mass spectrum.
}

\source{
The raw data were downloaded from this link : \url{https://zenodo.org/record/163517#.YIkWiNZuJCp}. The dataset focuses only on mass spectra from \emph{Citrobacter}.
}

\usage{
data("CitrobacterRKIspectra", package = "MSclassifR")

#####
#Plotting the first mass spectrum
#library("MSclassifR")
#PlotSpectra(SpectralData=CitrobacterRKIspectra[[1]],absx = "ALL", Peaks = NULL, 
#            Peaks2 = NULL, col_spec = 1, col_peak = 2, shape_peak = 3,
#            col_peak2 = 2, shape_peak2 = 2)
}

\description{
Mass spectra of the \code{\link{CitrobacterRKIspectra}} dataset.
}

\details{
The Robert Koch-Institute (RKI) database of microbial MALDI-TOF mass spectra contains raw mass spectra. Only mass spectra of the \emph{Citrobacter} bacterial species were collected.
}

\references{
Lasch, Peter, Stammler, Maren, & Schneider, Andy. (2018). Version 3 (20181130) of the MALDI-TOF Mass Spectrometry Database for Identification and Classification of Highly Pathogenic Microorganisms from the Robert Koch-Institute (RKI) [Data set]. Zenodo.\doi{http://doi.org/10.5281/zenodo.163517}
}

\keyword{datasets}
