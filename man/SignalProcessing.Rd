\name{SignalProcessing}

\alias{SignalProcessing}

\title{Function performing post acquisition signal processing}

\description{
This function performs post acquisition signal processing for \code{list} of \code{MassSpectrum} objects using commonly used methods : transform intensities ("sqrt"), smoothing ("Wavelet"), remove baseline ("SNIP"), calibrate intensities ("TIC") and align spectra. Methods used are selected from the \code{MALDIquant} and \code{MALDIrppa} R packages.
}

\usage{

SignalProcessing(x,
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
                 minFrequency= 0.5,
                 binPeaks_method = "strict",
                 keepReferenceSpectra = FALSE,
                 ...)
}

\arguments{
  \item{x}{a \code{list} of \code{MassSpectrum} objects (see \code{MALDIquant} R package).}

  \item{transformIntensity_method}{a \code{character} indicating the method used to transform intensities: \code{"sqrt"} by default. This function can be replaced by another mathematical function such as \code{"log"}.}

  \item{smoothing_method}{a \code{character} indicating the smoothing methods used. By default, it performs undecimated \code{Wavelet} transform (UDWT) for \code{list} of \code{MassSpectrum} objects. This Smoothing method can be remplaced by \code{"SavitzkyGolay"} or \code{"MovingAverage"}. See \code{wavSmoothing} in the \code{MALDIrppa} R package for details.}

  \item{removeBaseline_method}{a \code{character} indicating the method used to remove baseline. It uses \code{"SNIP"} method for \code{list} of \code{MassSpectrum} objects. This baseline estimation method can be remplaced \code{"TopHat"}, \code{"ConvexHull"} or \code{"median"}. See \code{removeBaseline-methods} of the \code{MALDIquant} R package for details.}

  \item{removeBaseline_iterations}{a \code{numeric} value indicting the number of iterations to remove baseline (by default = {\code{25}}). See \code{removeBaseline-methods} of the \code{MALDIquant} R package for details.}

  \item{calibrateIntensity_method}{a \code{character} indicating the intensities calibration method used (\code{"TIC"} method by default).  This calibration method can be remplaced by \code{"PQN"} or {"median"}.See \code{calibrateIntensity-methods} of the \code{MALDIquant} R package for details.}

   \item{alignSpectra_NoiseMethod}{a \code{character} indicating the noise estimation method. It uses \code{"MAD"} method for \code{list} of \code{MassSpectrum} objects. This noise estimation method estimation method can be remplaced \code{"SuperSmoother"}. See \code{estimateNoise-methods} of the \code{MALDIquant} R package for details.}
  
  \item{alignSpectra_method}{a \code{character} indicating the warping method. It uses \code{"lowess"} method for \code{list} of \code{MassSpectrum} objects. This warping method method can be remplaced \code{"linear"}, \code{"quadratic"} or \code{"cubic"} . See \code{determineWarpingFunctions} of the \code{MALDIquant} R package for details.}
  
  \item{alignSpectra_halfWs}{a \code{numeric} value half window size to detect peaks (by default = {\code{11}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{alignSpectra_SN}{a \code{numeric} value indicating the signal-to-noise ratio used to detect peaks (by default = {\code{3}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}
  
  \item{tolerance_align}{a \code{numeric} value indicating a maximal relative deviation of a peak position (mass) to be considered as identical in ppm (by default = \code{0.002}). See \code{determineWarpingFunctions} of the \code{MALDIquant} R package for details.}

  \item{referenceSpectra}{a \code{MassPeaks} reference spectrum for alignment of the sample(s). If this reference spectrum is not provided, a reference spectrum is created using the \code{minFrequency} and \code{binPeaks_method} arguments. See \code{referencePeaks} of the \code{MALDIquant} R package for details.}
  
  \item{minFrequency}{a \code{numeric} minimum frequency for each peak over all analyzed spectra (by default = \code{0.5}) for the creation of the reference spectrum. See \code{referencePeaks} of the \code{MALDIquant} R package for details.}
  
  \item{binPeaks_method}{a \code{character} indicating the method used to equalize masses for similar peaks for the creation of the reference spectrum. The \code{"strict"} method is used by default corresponding to a unique peak per bin from the same sample. This method can be remplaced by \code{"relaxed"} corresponding to multiple peaks per bin from the same sample. See \code{binPeaks} of the \code{MALDIquant} R package for more details.}
  
  \item{keepReferenceSpectra}{a \code{logical} value indicating if the created reference spectrum is returned by the function (\code{FALSE}, default).}

  \item{...}{other arguments from \code{MALDIrppa} packages for the \code{wavSmoothing} function such as \code{n.levels} (corresponding to the depth of the decomposiion for the wavelet function). See \code{wavSmoothing} of the \code{MALDIrppa} R package for details.}

}

\value{A list of modified \code{MassSpectrum} objects (see \code{MALDIquant} R package) according to chosen arguments. If the argument \code{referenceSpectra} is not completed and the argument \code{keepReferenceSpectra} is \code{TRUE}, a list containing the \code{MassSpectrum} objects modified named \code{"spectra"} and the created reference spectrum named \code{"RefS"} is returned.}

\references{

Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. Bioinformatics. 2012 Sep 1;28(17):2270-1. \doi{10.1093/bioinformatics/bts447}. Epub 2012 Jul 12. PMID: 22796955.

Javier Palarea-Albaladejo, Kevin Mclean, Frank Wright, David G E Smith, MALDIrppa: quality control and robust analysis for mass spectrometry data, Bioinformatics, Volume 34, Issue 3, 01 February 2018, Pages 522 - 523, \doi{https://doi.org/10.1093/bioinformatics/btx628}

}


\details{

The \code{SignalProcessing} function provides an analysis pipeline for \code{MassSpectrum} objects including intensity transformation, smoothing, removing baseline.

The \code{Wavelet} method relies on the \code{wavShrink} function of the \code{wmtsa} package and its dependencies (now archived by CRAN). The original C code by William Constantine and Keith L. Davidson, in turn including copyrighted routines by Insightful Corp., has been revised and included into \code{MALDIrppa} for the method to work.

All the methods used for \code{SignalProcessing} functions are selected from \code{MALDIquant} and \code{MALDIrppa} packages.
}


\examples{

library("MALDIquant")
library("MSclassifR")

## Load mass spectra
data("CitrobacterRKIspectra", package = "MSclassifR")

# plot first unprocessed mass spectrum
PlotSpectra(SpectralData=CitrobacterRKIspectra[[1]], col_spec="blue")

## spectral treatment
spectra <- SignalProcessing(CitrobacterRKIspectra,
                            transformIntensity_method = "sqrt",
                            smoothing_method = "Wavelet",
                            removeBaseline_method = "SNIP",
                            removeBaseline_iterations = 25,
                            calibrateIntensity_method = "TIC",
                            alignSpectra_Method = "MAD",
                            alignSpectra_halfWs = 11,
                            alignSpectra_SN = 3,
                            tolerance_align = 0.002)


# plot first processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]], col_spec="blue")

}

\keyword{Preprocessing}




