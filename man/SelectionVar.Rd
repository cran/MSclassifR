\name{SelectionVar}

\alias{SelectionVar}

\title{Variable selection using random forests, logistic regression methods or sparse partial least squares discriminant analysis (sPLS-DA).}

\description{
This function performs variable selection (i.e. selection of discriminant mass-over-charge values) using either recursive feature elimination (RFE) algorithm with Random Forest, or logistic regression model, or sparse partial least squares discriminant analysis (sPLS-DA).
}

\usage{
SelectionVar(X,
             Y,
             MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA"),
             MethodValidation = c("cv", "repeatedcv", "LOOCV"),
             PreProcessing = c("center","scale","nzv","corr"),
             NumberCV = NULL,
             RepeatsCV = NULL,
             Sizes,
             Ntree = 1000,
             threshold = 0.01,
             ncomp.max = 10
             )

}

\arguments{
  \item{X}{a numeric \code{matrix} corresponding to a library of mass spectra. Each row of \code{X} is the intensities of a mass spectrum measured on mass-over-charge values.}

 \item{Y}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}

  \item{MethodSelection}{a \code{charater} indicating the method used for variables selection. Four methods are avaible: with recursive feature elimination (RFE) and random forest (\code{"RFERF"}); logistic regression method (\code{"RFEGlmnet"}); method with random forest (\code{"VSURF"}) and sparse partial least squares discriminant analysis (\code{"sPLSDA"}).}

  \item{MethodValidation}{a \code{charater} indicating the resampling method:\code{"cv"} for cross-validation; \code{"repeatedcv"} for repeated cross-validation; and \code{"LOOCV"} for leave-one-out cross-validation.}

  \item{NumberCV}{a \code{numeric} value indicating K-folds for cross-validation. Don't used for \code{"VSURF"} method.}

   \item{RepeatsCV}{a \code{numeric} value indication the number of repeat(s) for K-folds for cross-validation or repeated cross-validation. Don't used for \code{"VSURF"} method.}

   \item{PreProcessing}{a \code{vector} indicating the method(s) used to pre-process the mass spectra in \code{X}: centering (\code{"center"}), scaling (\code{"scale"}), eliminating near zero variance predictors (\code{"nzv"}), or correlated predictors (\code{"corr"}).}

  \item{Sizes}{a numeric \code{vector} indicating the number of variables to select. Don't used for \code{"VSURF"} method.}

  \item{Ntree}{a \code{numeric} value indicating the number of trees in each forest, only used if \code{MethodSelection} = \code{"VSURF"} (1000 by default).}

  \item{ncomp.max}{a \code{positive Integer} indicating the maximum number of components included in the sPLS-DA model (10 by default).}

  \item{threshold}{a \code{positive Integer} corresponding to a threshold used for optimal selection of the number of components included in the sPLS-DA model (0.01 by default).}

}

\value{A list composed of:

\item{sel_moz}{a \code{vector} with discriminant mass-over-chage values.}

And of the results of the \code{rfe} function of the \code{caret} R package (methods \code{RFERF} and \code{RFEGlmnet}), or of the \code{VSURF} function of the \code{VSURF} R package (method \code{VSURF}). 

For the \code{sPLSDA} method,  it also returns the following items:
  \item{Raw_data}{a horizontal bar plot and  containing the contribution of features on each component.}

  \item{selected_variables}{\code{data frame} with uniques features (selected variables to keep and containing the contribution of features in order to class samples).See \code{plotLoadings} in the \code{mixOmics} R package for details.}
}


\details{
See  \code{rfe} in the \code{caret} R package, \code{VSURF} in the \code{VSURF} R package and \code{splsda} in the \code{mixOmics} R package for details.
}

\references{

Kuhn, Max. (2012). The caret Package. Journal of Statistical Software. 28.

Genuer, Robin, Jean-Michel Poggi and Christine Tuleau-Malot. VSURF : An R Package for Variable Selection Using Random Forests. R J. 7 (2015): 19.

Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22.

Kim-Anh Le Cao, Florian Rohart, Ignacio Gonzalez, Sebastien Dejean with key contributors Benoit Gautier, Francois, Bartolo, contributions from Pierre Monget, Jeff Coquery, FangZou Yao and Benoit Liquet. (2016). mixOmics: Omics. Data Integration Project. R package version 6.1.1. \url{https://CRAN.R-project.org/package=mixOmics}

}

\seealso{
\bold{Vignettes MSclassifR :}
\url{https://agodmer.github.io/MSclassifR_examples/Vignettes/Vignettemsclassifr_Ecrobia.html}
\url{https://agodmer.github.io/MSclassifR_examples/Vignettes/Vignettemsclassifr_Klebsiella.html}
}

\examples{

\donttest{

library("MSclassifR")
library("MALDIquant")

####################################
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


####################################
## 2. Perform variables selection using SelectionVar with RFE and random forest 
## (with 5 to 10 variables)
a <- SelectionVar(X,
                  Y,
                  MethodSelection = c("RFERF"),
                  MethodValidation = c("cv"),
                  PreProcessing = c("center","scale","nzv","corr"),
                  NumberCV = 2,
                  Sizes = c(5:10))
                  
# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
Peaks2=a$sel_moz,col_spec="blue",col_peak="black")

## 3. Perform variables selection unsing SelectionVar with RFE and logistic 
## regression (with 5 to 10 variables)
## It is recommended to have a large enough data set to use this method
UpFeatures <-caret::upSample(X,Y, list = TRUE)
b <- SelectionVar(UpFeatures$x,
                  UpFeatures$y,
                  MethodSelection = c("RFEGlmnet"),
                  MethodValidation = c("cv"),
                  PreProcessing = c("center","scale","nzv","corr"),
                  NumberCV = 2,
                  Sizes = c(5:10))

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
Peaks2=b$sel_moz,col_spec="blue",col_peak="black")

## 4. Perform variables selection unsing sPLDA method (with 5 to 10 variables per components)
#c <- SelectionVar(X,
#                  Y,
#                  MethodSelection = c("sPLSDA"),
#                  MethodValidation = c("LOOCV"),
#                  PreProcessing = c("scale","nzv"),
#                  Sizes = c(5:10))
                  
# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
#PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
#Peaks2=c$sel_moz,col_spec="blue",col_peak="black")

## 5. Perform variables selection using SelectionVar with RFE and logistic 
## regression (with 5 to 10 variables per components)
  ## This function can last a few minutes

  d <- SelectionVar(X, Y, MethodSelection = c("VSURF"))
  summary(d$result)

}

}
