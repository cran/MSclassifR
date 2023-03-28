\name{SelectionVarStat}

\alias{SelectionVarStat}

\title{Variable selection using statistical tests. Estimating the number of discriminant features (mass-over-chage values).}

\description{This function performs statistical tests for each mass-over-chage value to determine which are discriminants between categories. Using the distribution of resulting p-values, it determines the expected number of discriminant features.}


\usage{

SelectionVarStat(X,
                 Y,
                 stat.test = "Limma",
                 pi0.method="abh",
                 fdr=0.05,
                 Sampling = c("no", "up","down", "smote"))

}

\arguments{
 \item{X}{a numeric \code{matrix} corresponding to a library of mass spectra. Rows of \code{X} are the intensities of a mass spectrum measured on mass-over-charge values. The columns are mass-over-charge values.}

 \item{Y}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}
 
 \item{stat.test}{a \code{character} among "anova", "kruskal", or "Limma" (default). It corresponds to the test used to know if the intensity measured at a mass-over-charge value is significantly different between categories. "anova" is for a classical ANOVA Fisher test, "kruskal" is for the Kruskal-Wallis test, "Limma" is for an ANOVA Fisher test using the \code{limma} R package.}
 
 \item{pi0.method}{a \code{character} among "abh", "st.spline", "st.boot", "langaas", "histo", "pounds", "jiang", "slim". It corresponds to statistical methods used to estimate the proportion of true null hypotheses among the set of tested mass-over-charge values. See the \code{estim.pi0} function of th R package \code{cp4p} for details.}
 
 \item{fdr}{a \code{numeric} corresponding to False Discovery Rate threshold used to determine the differential mass-over-charge values. 0.05 by default.}
 
 \item{Sampling}{a \code{character} indicating an optional subsampling method to handle imbalanced datasets: subsampling methods are either \code{"no"} (no subsampling), \code{"up"}, \code{"down"} or \code{"smote"}. \code{"no"} by default.}
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
peaks <- MSclassifR::PeakDetection(x = spectra, averageMassSpec=FALSE)
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
## 2. Perform variables selection using SelectionVar with RFE and random forest 
# with 5 to 10 variables, 
# up sampling method and trained with the Kappa coefficient metric
a <- SelectionVar(X,
                  Y,
                  MethodSelection = c("RFERF"),
                  MethodValidation = c("cv"),
                  PreProcessing = c("center","scale","nzv","corr"),
                  NumberCV = 2,
                  Metric = "Kappa",
                  Sizes = c(5:10),
                  Sampling = "up")

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
            Peaks2=a$sel_moz,col_spec="blue",col_peak="black")

###############################################################################
## 3. Perform variables selection using SelectionVar with VSURF 
# This function can last a few minutes
d <- SelectionVar(X, Y, MethodSelection = c("VSURF"))
summary(d$result)

###############################################################################
## 4. Perform variables selection using SelectionVar with "mda" or "cvp"
# option 1: Using mean decrease in accuracy  
# with no sampling method
e <- SelectionVar(X,Y,MethodSelection="mda",Ntree=10*ncol(X)) 

# Estimation of the number of peaks to discriminate species
e$nb_to_sel

# Discriminant mass-over-charge values 
e$sel_moz

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
            Peaks2=e$sel_moz,col_spec="blue",col_peak="black")

# option 2: Using cross-validated permutation variable importance measures (more "time-consuming")       
# with no sampling method
f <- SelectionVar(X,Y,MethodSelection="cvp",NumberCV=2,ncores=2,Ntree=1000)

# Estimation of the number of peaks to discriminate species
f$nb_to_sel

# Discriminant mass-over-charge values 
f$sel_moz

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
            Peaks2=f$sel_moz,col_spec="blue",col_peak="black")

# Mass-over charge values found with both methods ("mda" and "cvp")
intersect(e$sel_moz,f$sel_moz)

}
}

\keyword{Feature selection}

