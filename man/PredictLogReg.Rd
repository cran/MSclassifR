\name{PredictLogReg}

\alias{PredictLogReg}

\title{Prediction of the category to which a mass spectrum belongs from a multinomial logistic regression model}

\description{
This function predicts the category (species, phenotypes...) to which a mass spectrum belongs from a set of shortlisted mass-over-charge values of interest and a short-listed multinomial logistic regression model (see \code{\link{LogReg}}).
}

\usage{
PredictLogReg(peaks,
               model,
               moz,
               tolerance = 6,
               toleranceStep = 2,
               normalizeFun = TRUE,
               noMatch=0,
               Reference = NULL)
}

\arguments{
  \item{peaks}{a list of \code{MassPeaks} objects (see \code{MALDIquant} R package).}

  \item{model}{a model or a list of models estimated from a set of shortlisted mass-over-charge values (output of the \code{\link{LogReg}} function).}

  \item{moz}{a \code{vector} with the set of shortlisted mass-over-charge values used to estimate the model \code{Model}.}

  \item{tolerance}{a \code{numeric} value of accepted tolerance to match peaks to the set of shortlisted mass-over-charge values. It is fixed to 6 Da by default.}
  
   \item{toleranceStep}{a \code{numeric} value added to the \code{tolerance} parameter to match peaks to the set of shortlisted mass-over-charge values. It is fixed to 2 Da by default.}

  \item{normalizeFun}{a \code{logical} value, if \code{TRUE} (default) the maximum intensity will be equal to 1, the other intensities will be expressed in ratio to this maximum.}

  \item{noMatch}{a \code{numeric} value used to replace intensity values if there is no match detected between peaks and the set of shortlisted mass-over-charge values \code{moz}. It is fixed to 0 by default.}

   \item{Reference}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}. \code{"NULL"} by default.}
 
}

\value{Returns a \code{dataframe} containing probabilities of membership by category for each mass spectrum in \code{peaks}. The method used is provided in the \code{method} column. The \code{comb_fisher} method is the result of the Fisher's method when merging probabilities of membership of used prediction models.The \code{max_vote} method is the result of the maximum voting from used prediction models.

If the \code{Reference} parameter is not null, the function returns:

  \item{Confusion.Matrix}{a \code{list} of confusion matrix (cross-tabulation with associated statitics) corresponding to the output of the \code{confusionMatrix} function of the \code{caret} R package.}

  \item{Gobal.stat}{a \code{data.frame} with three columns corresponding to the value (\code{value} column) of a statistic parameter (\code{Statistic.parameter} column) from a method used (\code{model} column) obtained with the \code{LogReg} function. See \code{LogReg} function for the Statistic.parameter column.}
  
  \item{Details.stat}{a \code{data.frame} with four columns corresponding to the same as \code{Gobal.stat} dataframe with the class concerned for estimated statistic parameter (class column). All statistic parameters are extracted from the output of the \code{confusionMatrix} function of the \code{caret} R package.}
  
  \item{Correct.ClassificationFreq}{a \code{data.frame} with predicted class (Prediction column) from a method (Model column) and the reference of the categories of each mass spectrum (Reference column). The \code{Freq} column indicates the number of times the category was correctly predicted by the method.}
  
  \item{Incorrect.ClassificationFreq}{a \code{data.frame} with predicted class (Prediction column) from a method (\code{Model} column) and the reference of the categories of each mass spectrum (\code{Reference} column). The \code{Freq} column indicates the number of times the category was not correctly predicted by the method. }

}

\details{The \code{PredictLogReg} function allows to predict the membership of a mass spectrum to a category from a multinomial logistic regression model. The mass spectrum from the \code{peaks} object will be matched to the discriminant mass-over-chage (m/z) values (\code{sel_moz} object from the \code{SelectionVar} function)  with a tolerance between two m/z defined by the \code{tolerance} parameter (by default this value is 6 Da). If a repetition of same m/z occurs in the selection, only the m/z that is closest in mass peaks (\code{moz}) is used. When no match, intensity values are replaced by the \code{noMatch} argument. If no m/z values from \code{peaks} object matched with the m/z in the object \code{moz}, the tolerance will be increased according to a numeric value defined in the \code{toleranceStep} parameter and a warning will be notified. Note that it is possible to not perform the \code{SelectionVar} function prior to the \code{PredictLogReg} function, and to replace the argument \code{moz} by all the m/z values present in a mass spectrum.}

\references{

Kuhn, M. (2008). Building predictive models in R using the caret package. Journal of statistical software, 28(1), 1-26.

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
spectra <- SignalProcessing(CitrobacterRKIspectra)
# detection of peaks in pre-processed mass spectra
peaks <- PeakDetection(x = spectra, labels = CitrobacterRKImetadata$Strain_name_spot)
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
## 2. Selection of discriminant mass-over-charge values using sPLS-DA
# with 5 to 10 variables, 
# up-sampling method and trained 
# with the Accuracy coefficient metric

a <- MSclassifR::SelectionVar(X,
                              Y,
                              MethodSelection = c("RFERF"),
                              MethodValidation = c("cv"),
                              PreProcessing = c("center","scale","nzv","corr"),
                              NumberCV = 2,
                              Metric = "Accuracy",
                              Sizes = c(2:5),
                              Sampling = "up")

sel_moz=a$sel_moz

###############################################################################
## 3. Perform LogReg from shortlisted discriminant mass-over-charge values

# linear multinomial regression 
# without sampling mehod and 
# trained with the Kappa coefficient metric

model_lm=MSclassifR::LogReg(X=X,
                            moz=sel_moz,
                            Y=factor(Y),
                            number=2,
                            repeats=2,
                            Metric = "Kappa")
# Estimated model:
model_lm

# nonlinear multinomial regression using neural networks 
# with up-sampling method and 
# trained with the Kappa coefficient metric

model_nn=MSclassifR::LogReg(X=X,
                            moz=sel_moz,
                            Y=factor(Y),
                            number=2,
                            repeats=2,
                            kind="nnet",
                            Metric = "Kappa",
                            Sampling = "up")
# Estimated model:
model_nn

# nonlinear multinomial regression using random forests 
# without down-sampling method and 
# trained with the Kappa coefficient metric

model_rf=MSclassifR::LogReg(X=X,
                            moz=sel_moz,
                            Y=factor(Y),
                            number=2,
                            repeats=2,
                            kind="rf",
                            Metric = "Kappa",
                            Sampling = "down")

# Estimated model:
model_rf

# nonlinear multinomial regression using xgboost 
# with down-sampling method and 
# trained with the Kappa coefficient metric

model_xgb=MSclassifR::LogReg(X=X,
                             moz=sel_moz,
                             Y=factor(Y),
                             number=2,
                             repeats=2,
                             kind="xgb",
                             Metric = "Kappa",
                             Sampling = "down")
# Estimated model:
model_xgb

# nonlinear multinomial regression using svm 
# with down-sampling method and 
# trained with the Kappa coefficient metric

model_svm=MSclassifR::LogReg(X=X,
                             moz=sel_moz,
                             Y=factor(Y),
                             number=2,
                             repeats=2,
                             kind="svm",
                             Metric = "Kappa",
                             Sampling = "down")
# Estimated model:
model_svm

# Of note, you can also load a model already saved 
# (see example in LogReg function) for the next step
###############################################################################
## 4. Probabilities of belonging to each category for the mass spectra 
## and associated statitics

# Collect all the estimated models in a list

Models <- list(model_lm$train_mod,
               model_nn$train_mod,
               model_rf$train_mod,
               model_xgb$train_mod,
               model_svm$train_mod)

# Predict classes of mass spectra with 6 Da of tolerance for matching peaks.
prob_cat=MSclassifR::PredictLogReg(peaks = peaks[c(1:5)],
                                   model = Models,
                                   moz = sel_moz,
                                   tolerance = 6,
                                   Reference = Y[c(1:5)])

}
}

\keyword{Prediction}