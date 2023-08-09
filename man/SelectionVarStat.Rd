\name{SelectionVarStat}

\alias{SelectionVarStat}

\title{Variable selection using multiple statistical tests.}

\description{
This function performs a statistical test for each mass-over-charge value to determine which are discriminants between categories. Using the distribution of resulting multiple p-values, it determines an expected number of discriminant features, and adjusted p-values that can be used to control a false discovery rate threshold.
}


\usage{

SelectionVarStat(X,
                 Y,
                 stat.test = "Limma",
                 pi0.method="abh",
                 fdr=0.05,
                 Sampling = c("no", "up","down", "smote"))

}

\arguments{
 \item{X}{a \code{numeric} matrix corresponding to a library of mass spectra. Rows of \code{X} are the intensities of a mass spectrum measured on mass-over-charge values. The columns are mass-over-charge values.}

 \item{Y}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}
 
 \item{stat.test}{a \code{character} among \code{"anova"}, \code{"kruskal"}, or \code{"Limma"} (default). It corresponds to the test used to know if the intensity measured at a mass-over-charge value is significantly different between categories. \code{"anova"} is for a classical ANOVA Fisher test, \code{"kruskal"} is for the Kruskal-Wallis test, \code{"Limma"} is for an ANOVA Fisher test using the \code{limma} R package.}
 
 \item{pi0.method}{a \code{character} among \code{"abh"}, \code{"st.spline"}, \code{"st.boot"}, \code{"langaas"}, \code{"histo"}, \code{"pounds"}, \code{"jiang"}, \code{"slim"}. It corresponds to statistical methods used to estimate the proportion of true null hypotheses among the set of tested mass-over-charge values. See the \code{estim.pi0} function of the R package \code{cp4p} for details.}
 
 \item{fdr}{a \code{numeric} value corresponding to False Discovery Rate threshold used to determine the differential mass-over-charge values. 0.05 by default.}
 
 \item{Sampling}{a \code{character} indicating an optional subsampling method to handle imbalanced datasets: subsampling methods are either \code{"no"} (no subsampling), \code{"up"}, \code{"down"} or \code{"smote"}. \code{"no"} by default.}
}

\value{A list composed of:

  \item{nb_to_sel}{a \code{numeric} value corresponding to an estimation of the optimal number of mass-over-charge values to discriminate between different groups.}
  
  \item{sel_moz}{a \code{vector} with selected discriminant mass-over-charge values.}
  
  \item{ap}{a \code{list} composed of \code{pi0} the proportion of non-discriminant mass-over-charge values, and \code{adjp} a matrix of raw p-values and corresponding ajusted p-values for all the mass-over-charge values that have been tested.}
}


\details{The \code{SelectionVarStat} function allows performing "quick" classification of mass-over-charge values. It tries to find all the mass-over-charge values (or the number of mass-over-charge values) that are discriminant between categories. This can conduct to select "correlated" mass-over-charge values (i.e. associated to intensities evolving similarly between categories). By default, multiple moderated t-tests using the \code{limma} R package (bayesian regularization of variances) are performed and the p-values are corrected using an adaptive Benjamini and Hochberg procedure to control the false discovery rate. Different ways to estimate the proportion of true null hypotheses (object \code{pi0} returned by the function - see the  \code{cp4p} R package for details) can be used for the adaptive Benjamini-Hochberg procedure ("\code{abh}" by defaut).}

\references{

Gianetto, Quentin & Combes, Florence & Ramus, Claire & Bruley, Christophe & Coute, Yohann & Burger, Thomas. (2015). Technical Brief Calibration Plot for Proteomics (CP4P): A graphical tool to visually check the assumptions underlying FDR control in quantitative experiments. Proteomics. 16. 10.1002/pmic.201500189. 

}

\examples{
\donttest{
library("MSclassifR")
library("MALDIquant")

###############################################################################
## 1. Pre-processing of mass spectra

# load mass spectra and their metadata
data("CitrobacterRKIspectra","CitrobacterRKImetadata", package = "MSclassifR")
# standard pre-processing of mass spectra
spectra <- MSclassifR::SignalProcessing(CitrobacterRKIspectra)
# detection of peaks in pre-processed mass spectra
peaks <- MSclassifR::PeakDetection(x = spectra, labels = CitrobacterRKImetadata$Strain_name_spot)
# matrix with intensities of peaks arranged in rows (each column is a mass-over-charge value)
IntMat <- MALDIquant::intensityMatrix(peaks)
rownames(IntMat) <- paste(CitrobacterRKImetadata$Strain_name_spot)
# remove missing values in the matrix
IntMat[is.na(IntMat)] <- 0
# normalize peaks according to the maximum intensity value for each mass spectrum
IntMat <- apply(IntMat,1,function(x) x/(max(x)))
# transpose the matrix for statistical analysis
X <- t(IntMat)
# define the known categories of mass spectra for the classification
Y <- factor(CitrobacterRKImetadata$Species)

###############################################################################
## 2. Estimate the optimal number of peaks to discriminate the different species

OptiPeaks <- SelectionVarStat(X,
                              Y,
                              stat.test = "Limma",
                              pi0.method="abh",
                              fdr=0.05,
                              Sampling="smote")
             
## Estimation of the optimal number of peaks to discriminate species (from the pi0 parameter)
OptiPeaks$nb_to_sel

## discriminant mass-over-chage values estimated using a 5 per cent false discovery rate
OptiPeaks$sel_moz

## p-values and adjusted p-values estimated for all the tested mass-over-charge values
OptiPeaks$ap$adjp

}
}

\keyword{Feature selection}

