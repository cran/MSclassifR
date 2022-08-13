\name{SelectionVarStat}

\alias{SelectionVarStat}

\title{Variable selection using statistical tests. Estimating the number of discriminant features (mass-over-chage values).}

\description{This function performs statistical tests for each mass-over-chage value to determine which are discriminants between categories. Using the distribution of resulting p-values, it determines the expected number of discriminant features.}


\usage{

SelectionVarStat(X,
                 Y,
                 stat.test = "Limma",
                 pi0.method="abh",
                 fdr=0.05)

}

\arguments{
 \item{X}{a numeric \code{matrix} corresponding to a library of mass spectra. Rows of \code{X} are the intensities of a mass spectrum measured on mass-over-charge values. The columns are mass-over-charge values.}

 \item{Y}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}
 
 \item{stat.test}{a \code{character} among "anova", "kruskal", or "Limma" (default). It corresponds to the test used to know if the intensity measured at a mass-over-charge value is significantly different between categories. "anova" is for a classical ANOVA Fisher test, "kruskal" is for the Kruskal-Wallis test, "Limma" is for an ANOVA Fisher test using the \code{limma} R package.}
 
 \item{pi0.method}{a \code{character} among "abh", "st.spline", "st.boot", "langaas", "histo", "pounds", "jiang", "slim". It corresponds to statistical methods used to estimate the proportion of true null hypotheses among the set of tested mass-over-charge values. See the \code{estim.pi0} function of th R package \code{cp4p} for details.}
 
 \item{fdr}{a \code{numeric} corresponding to False Discovery Rate threshold used to determine the differential mass-over-charge values. 0.05 by default.}
}

\value{A list composed of:

  \item{nb_to_sel}{a \code{numeric} value corresponding to an estimated number of mass-over-chage values where the intensities are significantly different between categories. It depends on the statistical methods used in \code{pi0.method}.}
  
  \item{NbEstimatedPeaks}{a \code{vector} with discriminant mass-over-chage values resulting to the FDR threshold applied on the set of tested mass-over-charge values.}
}


\details{The \code{SelectionVarStat} function allows performing "quick" classification of mass-over-charge values. It tries to find all the mass-over-charge values (or the number of mass-over-charge values) that are discriminant between categories. This can conduct to select "correlated" (i.e. associated to intensities evolving similarly between categories) mass-over-charge values.}

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
## 2. Perform the SelectionVarStat function
                
OptiPeaks <- SelectionVarStat(X,
                              Y,
                              stat.test = ("Limma"),
                              pi0.method="abh",
                              fdr=0.05)
             
## Estimation of the number of peaks to discriminate species
OptiPeaks$nb_to_sel

## Discriminant mass-over-chage values 
OptiPeaks$sel_moz

}
}



