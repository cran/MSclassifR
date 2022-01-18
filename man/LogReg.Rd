\name{LogReg}

\alias{LogReg}

\title{Estimation of a multinomial logistic regression to predict the category to which a mass spectrum belongs}

\description{
This function estimates a multinomial logistic regression using cross-validation to predict the category (species, phenotypes...) to which a mass spectrum belongs from a set of shortlisted mass-over-charge values corresponding to discriminant peaks. Two main kinds of models can be estimated: linear or nonlinear (with neural networks, random forests, support vector machines with linear kernel, or eXtreme Gradient Boosting). Hyperparameters are randomly searched, except for the eXtreme Gradient Boosting where a grid search is performed.
}

\usage{
LogReg(X, moz, Y, number = 2, repeats = 2, kind="linear")
}

\arguments{
  \item{X}{\code{matrix} corresponding to a library of mass spectra. Each row of \code{X} is the intensities of a mass spectrum measured on the \code{moz} values.}

  \item{moz}{\code{vector} with shortlisted mass-over-charge values.}

  \item{Y}{\code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}

  \item{number}{\code{integer} corresponding to the number of folds or number of resampling iterations. See arguments of the \code{trainControl} function of the \code{caret} R package.}

  \item{repeats}{\code{integer} corresponding to the number of complete sets of folds to compute. See \code{trainControl} function of the \code{caret} R package for more details.}

  \item{kind}{If \code{kind="nnet"}, then a nonlinear multinomial logistic regression is estimated via neural networks. If \code{kind="rf"}, then it is estimated via random forests. If \code{kind="svm"}, then it is estimated via support vector machines with linear kernel. If \code{kind="xgb"}, then it is estimated via eXtreme gradient boosting. Else a linear multinomial logistic regression is performed (by default).}
}

\value{Returns a \code{list} with four items:

  \item{train_mod}{a \code{list} corresponding to the output of the train function of the \code{caret} R package containing the multinomial regression model estimated using repeated cross-validation.}

  \item{conf_mat}{a confusion matrix containing percentages classes of predicted categories in function of actual categories, resulting from repeated cross-validation.}

  \item{stats_global}{a \code{data frame} containing the mean and standard deviation values of the "Accuracy"" and "Kappa" parameters computed for each cross-validation.}

  \item{boxplot}{a \code{ggplot} object (see \code{ggplot2} R package). This is a graphical representation of the "Accuracy" and "Kappa"" parameters of \code{stats_global} using boxplots.}

}

\details{This function estimates a model from a library of mass spectra for which we already know the category to which they belong (ex.: species, etc). This model can next be used to predict the category of a new coming spectrum for which the category is unknown (see \code{\link{Predict_LogReg}}).
The estimation is performed using the \code{train} function of the \code{caret} R package. For each kind of model, random parameters are tested to find a model with the best Accuracy.
}

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

#linear multinomial regression
model_lm=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2)
#Estimated model:
model_lm

#nonlinear multinomial regression using neural networks
model_nn=MSclassifR::LogReg(X=X, moz=sel_moz, Y=factor(Y), number=2, repeats=2, kind="nnet")
#Estimated model:
model_nn

#nonlinear multinomial regression using random forests
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

##########
#Of note, step 3 can be performed several times to find optimal models 
#because of random hyperparameter search

#######################################
## 4. Select best models in term of average accuracy and saving it for reuse

acc_model=c(model_lm$stats_global[1,1],model_nn$stats_global[1,1],
model_rf$stats_global[1,1],model_xgb$stats_global[1,1],model_svm$stats_global[1,1])
names(acc_model)=c("lm","nn","rf","xgb","svm")
#Best models in term of accuracy
acc_model[which(acc_model==max(acc_model))]

#save best models for reuse
#models=list(model_lm$train_mod,model_nn$train_mod,model_rf$train_mod,
#model_xgb$train_mod,model_svm$train_mod)
#models_best=models[which(acc_model==max(acc_model))]
#for (i in 1:length(models_best)){
#save(models_best[[i]], file = paste0("model_best_",i,".rda",collapse="")
#}

#load a saved model
#load("model_best_1.rda")

}

}
