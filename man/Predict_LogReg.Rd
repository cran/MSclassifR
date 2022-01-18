\name{Predict_LogReg}

\alias{Predict_LogReg}

\title{Prediction of the category to which a mass spectrum belongs from a multinomial logistic regression model}

\description{
This function predicts the category (species, phenotypes...) to which a mass spectrum belongs from a set of shortlisted mass-over-charge values of interest and a short-listed multinomial logistic regression model (see \code{\link{LogReg}}).
}

\usage{
Predict_LogReg(peaks,model,moz,tolerance=6,normalizeFun=TRUE,noMatch=0)
}

\arguments{
  \item{peaks}{a list of \code{MassPeaks} objects (see \code{MALDIquant} R package).}

  \item{model}{a model or a list of models estimated from a set of shortlisted mass-over-charge values (output of the \code{\link{LogReg}} function).}

  \item{moz}{a \code{vector} with the set of shortlisted mass-over-charge values used to estimate the model \code{Model}.}

  \item{tolerance}{a \code{numeric} value of accepted tolerance to match peaks to the set of shortlisted mass-over-charge values. It is fixed to 6 m/z by default.}

  \item{normalizeFun}{a \code{logical} value, if \code{TRUE} (default) the maximum intensity will be equal to 1, the other intensities will be expressed in ratio to this maximum.}

  \item{noMatch}{a \code{numeric} value used to replace intensity values if there is no match detected between peaks and the set of shortlisted mass-over-charge values \code{moz}. It is fixed to 0 by default.}
}

\value{Returns a \code{dataframe} containing probabilities of membership by category for each mass spectrum in \code{peaks}. The method used is provided in the \code{method} column. The \code{comb_fisher} method is the result of the Fisher's method when merging probabilities of membership of used prediction models.The \code{max_vote} method is the result of the maximum voting from used prediction models.}

\references{
Kuhn, M. (2008). Building predictive models in R using the caret package. Journal of statistical software, 28(1), 1-26.
}

\examples{

\donttest{

library("MSclassifR")
library("MALDIquant")
library("mixOmics")

####################################
## 1. Pre-processing of mass spectra

# load mass spectra and their metadata
data("CitrobacterRKIspectra","CitrobacterRKImetadata", package = "MSclassifR")
# standard pre-processing of mass spectra
spectra <- SignalProcessing(CitrobacterRKIspectra)
# detection of peaks in pre-processed mass spectra
peaks <- PeakDetection(x = spectra, labels = CitrobacterRKImetadata$Number_strain)
# matrix with intensities of peaks arranged in rows (each column is a mass-over-charge value)
IntMat <- MALDIquant::intensityMatrix(peaks)
rownames(IntMat) <- paste(CitrobacterRKImetadata$Strain_name_spot)
# replace missing values with 0 in the matrix
IntMat[is.na(IntMat)] <- 0
# normalize peaks according to the maximum intensity value for each mass spectrum
IntMat <- apply(IntMat,1,function(x) x/(max(x)))
# transpose the matrix for statistical analysis
X <- t(IntMat)
# define the known categories of mass spectra for the classification
Y <- factor(CitrobacterRKImetadata$Species)

#####################################################################
## 2. Selection of discriminant mass-over-charge values using sPLS-DA

c <- SelectionVar(X,
                  Y,
                  MethodSelection = c("sPLSDA"),
                  MethodValidation = c("LOOCV"),
                  PreProcessing = c("scale","nzv"),
                  NumberCV = 2,
                  Sizes = c(5:20))
sel_moz=c$sel_moz

##########################################################################
## 3. Perform LogReg from shortlisted discriminant mass-over-charge values 
# and predict category of a mass spectrum

# Linear multinomial regression
model_lm=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2)

#nonlinear multinomial regression using neural networks
model_nn=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2, kind="nnet")
#Estimated model:
model_nn

#nonlinear multinomial regression using neural networks
model_rf=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2, kind="rf")
#Estimated model:
model_rf

#nonlinear multinomial regression using xgboost
model_xgb=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2, kind="xgb")
#Estimated model:
model_xgb

#nonlinear multinomial regression using svm
model_svm=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2, kind="svm")
#Estimated model:
model_svm

#Of note, you can also load a model already saved (see example in LogReg function)

## Probabilities of belonging to each category for the mass spectra
prob_cat=MSclassifR::Predict_LogReg(peaks = peaks[c(1:3)],
model = list(model_lm$train_mod,model_nn$train_mod,model_rf$train_mod,
model_xgb$train_mod,model_svm$train_mod),
moz = sel_moz)

}

}

