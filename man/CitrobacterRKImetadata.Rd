\docType{data}

\name{CitrobacterRKImetadata}

\alias{CitrobacterRKImetadata}

\title{Metadata of mass spectra corresponding to the bacterial species \emph{Citrobacter} sp. from  The Robert Koch-Institute (RKI) database of microbial MALDI-TOF mass spectra}

\format{
A data frame with 14 rows (each corresponding to a mass spectrum), and five columns that contain (in order): the strain name, the species name, the spot, a sample number and the name of the strain associated with the spot.
}

\source{
The raw data were downloaded from this link : \url{https://zenodo.org/record/163517#.YIkWiNZuJCp}. The dataset focuses only on mass spectra from \emph{Citrobacter}.
}

\usage{
data("CitrobacterRKImetadata", package = "MSclassifR")
}

\description{
Metadada of the \code{\link{CitrobacterRKIspectra}} list of mass spectra.
}

\details{
The Robert Koch-Institute (RKI) database of microbial MALDI-TOF mass spectra contains raw mass spectra. Only mass spectra of the \emph{Citrobacter} bacterial species were collected. Metadata were manually reported from raw data.
}

\references{
Lasch, Peter, Stammler, Maren, & Schneider, Andy. (2018). Version 3 (20181130) of the MALDI-TOF Mass Spectrometry Database for Identification and Classification of Highly Pathogenic Microorganisms from the Robert Koch-Institute (RKI) [Data set]. Zenodo.\doi{http://doi.org/10.5281/zenodo.163517}
}

\keyword{Dataset}
