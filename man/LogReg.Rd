\name{LogReg}

\alias{LogReg}

\title{Estimation of a multinomial regression to predict the category to which a mass spectrum belongs}

\description{
This function estimates a multinomial regression using cross-validation to predict the category (species, phenotypes...) to which a mass spectrum belongs from a set of shortlisted mass-to-charge values corresponding to discriminant peaks. Two main kinds of models can be estimated: linear or nonlinear (with neural networks, random forests, support vector machines with linear kernel, or eXtreme Gradient Boosting). Hyperparameters are randomly searched, except for the eXtreme Gradient Boosting where a grid search is performed.
}

\usage{
LogReg(X,
       moz,
       Y,
       number = 2,
       repeats = 2,
       Metric = c("Kappa", "Accuracy", "F1", "AdjRankIndex", "MatthewsCorrelation"),
       kind="linear",
       Sampling = c("no", "up", "down", "smote"))
}

\arguments{
  \item{X}{\code{matrix} corresponding to a library of mass spectra. Each row of \code{X} is the intensities of a mass spectrum measured on the \code{moz} values. The columns should be represented by mass-to-charge values.}

  \item{moz}{\code{vector} with shortlisted mass-to-charge values.}

  \item{Y}{\code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}

  \item{number}{\code{integer} corresponding to the number of folds or number of resampling iterations. See arguments of the \code{trainControl} function of the \code{caret} R package.}

  \item{Metric}{a \code{character} indicating metric to select the optimal model. The possibles metrics are the \code{"Kappa"} coefficient,\code{"Accuracy"}, the \code{"F1"} score, \code{"AdjRankIndex"} for the Adjusted Rand Index or \code{"MatthewsCorrelation"} for the Matthews Correlation Coefficient.}

  \item{repeats}{\code{integer} corresponding to the number of complete sets of folds to compute. See \code{trainControl} function of the \code{caret} R package for more details.}

  \item{kind}{If \code{kind="nnet"}, then a nonlinear multinomial logistic regression is estimated via neural networks. If \code{kind="rf"}, then it is estimated via random forests. If \code{kind="svm"}, then it is estimated via support vector machines with linear kernel. If \code{kind="xgb"}, then it is estimated via eXtreme gradient boosting. Else a linear multinomial logistic regression is performed (by default).}

  \item{Sampling}{a \code{character} indicating an optional subsampling method to handle imbalanced datasets: subsampling methods are either \code{"no"} (no subsampling), \code{"up"}, \code{"down"} or \code{"smote"}. \code{"no"} by default.}

}

\value{Returns a \code{list} with four items:

  \item{train_mod}{a \code{list} corresponding to the output of the train function of the \code{caret} R package containing the multinomial regression model estimated using repeated cross-validation.}

  \item{conf_mat}{a confusion matrix containing percentages classes of predicted categories in function of actual categories, resulting from repeated cross-validation.}

  \item{stats_global}{a \code{data frame} containing the mean and standard deviation values of the "Accuracy"" and "Kappa" parameters computed for each cross-validation.}

  \item{boxplot}{a \code{ggplot} object (see \code{ggplot2} R package). This is a graphical representation of the \code{Metric} parameters of \code{stats_global} using boxplots.}

}

\details{This function estimates a model from a library of mass spectra for which we already know the category to which they belong (ex.: species, etc). This model can next be used to predict the category of a new coming spectrum for which the category is unknown (see \code{\link{PredictLogReg}}).

The estimation is performed using the \code{train} function of the \code{caret} R package. For each kind of model, random parameters are tested to find a model according to the best \code{metric}. The formulas for the \code{metric} are as follows:

\deqn{Accuracy = Number Of Correct Predictions/Total Number Of Predictions}
\deqn{Kappa coefficient = (Observed Agreement-Chance Agreement)/(1-Chance Agreement)}
\deqn{F1 = True Positive/(True Positive + 1/2 (False Positive + False Negative))}

The adjusted Rand index (\code{"AdjRankIndex"}) is defined as the corrected-for-chance version of the Rand index which allows comparing two groups (see \code{mclust} package and \code{adjustedRandIndex()} function for more details). The Matthews correlation coefficient (\code{"MatthewsCorrelation"}) is estimated using \code{mcc} function in the \code{mltools} R package.

The \code{Sampling} methods available for imbalanced data are: \code{"up"} to the up-sampling method which consists of random sampling (with replacement) so that the minority class is the same size as the majority class; \code{"down"} to the down-sampling method which consists of random sampling (without replacement) of the majority class so that their class frequencies match the minority class; \code{"smote"} to the Synthetic Minority Over sampling Technique (SMOTE) algorithm for data augmentation which consists of creating new data from minority class using the K Nearest Neighbor algorithm.

Godmer et al. (2025) presents a comparison of different pipelines using LogReg that can help you to optimize your workflow. For a comprehensive guide, additional applications, and detailed examples of using this package, please visit our GitHub repository: \href{https://github.com/agodmer/MSclassifR_examples}{here}.
}

\references{
Kuhn, M. (2008). Building predictive models in R using the caret package. Journal of statistical software, 28(1), 1-26.

L. Hubert and P. Arabie (1985) Comparing Partitions, Journal of the Classification, 2, pp. 193-218.

Scrucca L, Fop M, Murphy TB, Raftery AE (2016). mclust 5: clustering, classification and density estimation using Gaussian finite mixture models. The R Journal.

Nitesh V. Chawla, Kevin W. Bowyer, Lawrence O. Hall, and W. Philip Kegelmeyer. 2002. SMOTE: synthetic minority over-sampling technique. J. Artif. Int. Res. 16, 1.

Matthews, B. W. (1975). "Comparison of the predicted and observed secondary structure of T4 phage lysozyme". Biochimica et Biophysica Acta (BBA) - Protein Structure. PMID 1180967.

Alexandre Godmer, Yahia Benzerara, Emmanuelle Varon, Nicolas Veziris, Karen Druart, Renaud Mozet, Mariette Matondo, Alexandra Aubry, Quentin Giai Gianetto, MSclassifR: An R package for supervised classification of mass spectra with machine learning methods, Expert Systems with Applications, Volume 294, 2025, 128796, ISSN 0957-4174, \doi{10.1016/j.eswa.2025.128796}.
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
peaks <- MSclassifR::PeakDetection(x = spectra, averageMassSpec=FALSE)
# matrix with intensities of peaks arranged in rows (each column is a mass-to-charge value)
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
## 2. Selection of discriminant mass-to-charge values using RFERF
# with 5 to 10 variables,
# up sampling method and
# trained with the Accuracy coefficient metric

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
## 3. Perform LogReg from shortlisted discriminant mass-to-charge values

# linear multinomial regression
# without sampling mehod
# and trained with the Kappa coefficient metric

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

##########
# Of note, step 3 can be performed several times
# to find optimal models
# because of random hyperparameter search

###############################################################################
## 4. Select best models in term of average Kappa and saving it for reuse

Kappa_model=c(model_lm$stats_global[1,2],model_nn$stats_global[1,2],
              model_rf$stats_global[1,2],model_xgb$stats_global[1,2],model_svm$stats_global[1,2])
names(Kappa_model)=c("lm","nn","rf","xgb","svm")
#Best models in term of accuracy
Kappa_model[which(Kappa_model==max(Kappa_model))]

#save best models for reuse
#models=list(model_lm$train_mod,model_nn$train_mod,model_rf$train_mod,
#model_xgb$train_mod,model_svm$train_mod)
#models_best=models[which(Kappa_model==max(Kappa_model))]
#for (i in 1:length(models_best)){
#save(models_best[[i]], file = paste0("model_best_",i,".rda",collapse="")
#}

#load a saved model
#load("model_best_1.rda")

###############################################################################
## 5. Try other metrics to select the best model

# linear multinomial regression
# with up-sampling method and
# trained with the Adjusted Rank index metric

model_lm=MSclassifR::LogReg(X=X,
                            moz=sel_moz,
                            Y=factor(Y),
                            number=2,
                            repeats=3,
                            Metric = "AdjRankIndex",
                            Sampling = "up")

}

}

\keyword{Estimation}
