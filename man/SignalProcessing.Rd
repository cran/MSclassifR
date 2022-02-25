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
                 alignSpectra_halfWs = 11,
                 alignSpectra_SN = 3,
                 tolerance_align = 0.002,
                 ...)
}

\arguments{
  \item{x}{a \code{list} of \code{MassSpectrum} objects (see \code{MALDIquant} R package).}

  \item{transformIntensity_method}{a \code{character} indicating the method used to transform intensities: \code{"sqrt"} by default.}

  \item{smoothing_method}{a \code{character} indicating the smoothing methods used. By default, it performs undecimated \code{Wavelet} transform (UDWT) for \code{list} of \code{MassSpectrum} objects. See \code{wavSmoothing} in the \code{MALDIrppa} R package for details.}

  \item{removeBaseline_method}{a \code{character} indicating the method used to remove baseline. It uses \code{"SNIP"} method for \code{list} of \code{MassSpectrum} objects. See \code{removeBaseline-methods} of the \code{MALDIquant} R package for details.}

  \item{removeBaseline_iterations}{a \code{numeric} value indicting the number of iterations to remove baseline (by default = {\code{25}}). See \code{removeBaseline-methods} of the \code{MALDIquant} R package for details.}

  \item{calibrateIntensity_method}{a \code{character} indicating the intensities calibration method used (\code{"TIC"} method by default). See \code{calibrateIntensity-methods} of the \code{MALDIquant} R package for details.}

  \item{alignSpectra_halfWs}{a \code{numeric} value half window size to detect peaks (by default = {\code{11}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{alignSpectra_SN}{a \code{numeric} value indicating the signal-to-noise ratio used to detect peaks (by default = {\code{3}}). See \code{detectPeaks-methods} of the \code{MALDIquant} R package for details.}

  \item{tolerance_align}{a \code{numeric} value indicating a maximal relative deviation of a peak position (mass) to be considered as identical (by default = \code{0.002}). See \code{determineWarpingFunctions} of the \code{MALDIquant} R package for details.}

  \item{...}{other arguments from \code{MALDIquant} and \code{MALDIrppa} packages.}

}

\value{A list of modified \code{MassSpectrum} objects (see \code{MALDIquant} R package) according to chosen arguments.}

\references{

Gibb S, Strimmer K. MALDIquant: a versatile R package for the analysis of mass spectrometry data. Bioinformatics. 2012 Sep 1;28(17):2270-1. \doi{10.1093/bioinformatics/bts447}. Epub 2012 Jul 12. PMID: 22796955.

Javier Palarea-Albaladejo, Kevin Mclean, Frank Wright, David G E Smith, MALDIrppa: quality control and robust analysis for mass spectrometry data, Bioinformatics, Volume 34, Issue 3, 01 February 2018, Pages 522 - 523, \doi{https://doi.org/10.1093/bioinformatics/btx628}

}


\details{
The \code{Wavelet} method relies on the \code{wavShrink} function of the \code{wmtsa} package and its dependencies (now archived by CRAN). The original C code by William Constantine and Keith L. Davidson, in turn including copyrighted routines by Insightful Corp., has been revised and included into \code{MALDIrppa} for the method to work.

All the methods used for \code{SpectralTreatment} functions are selected from \code{MALDIquant} and \code{MALDIrppa} packages.
}

\seealso{
\bold{Vignettes MSclassifR :}
\url{https://agodmer.github.io/MSclassifR_examples/Vignettes/Vignettemsclassifr_Ecrobia.html}
\url{https://agodmer.github.io/MSclassifR_examples/Vignettes/Vignettemsclassifr_Klebsiella.html}
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
                            alignSpectra_halfWs = 11,
                            alignSpectra_SN = 3,
                            tolerance_align = 0.002)


# plot first processed mass spectrum
PlotSpectra(SpectralData=spectra[[1]], col_spec="blue")

}



