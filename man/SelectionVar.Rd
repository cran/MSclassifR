\name{SelectionVar}

\alias{SelectionVar}

\title{Variable selection using random forests, logistic regression methods or sparse partial least squares discriminant analysis (sPLS-DA).}

\description{
This function performs variable selection (i.e. selection of discriminant mass-over-charge values) using either recursive feature elimination (RFE) algorithm with Random Forest, or logistic regression model, or sparse partial least squares discriminant analysis (sPLS-DA) or methods based on the distribution of variable importances of random forests.
}

\usage{
SelectionVar(X,
             Y,
             MethodSelection = c("RFERF", "RFEGlmnet", "VSURF", "sPLSDA", "mda", "cvp"), 
             MethodValidation = c("cv", "repeatedcv", "LOOCV"),
             PreProcessing = c("center", "scale", "nzv", "corr"), 
             Metric = c("Kappa", "Accuracy"), Sampling = c("no", "up","down", "smote"),
             NumberCV = NULL,
             RepeatsCV = NULL, 
             Sizes,
             Ntree = 1000,
             ncores = 2,
             threshold = 0.01,
             ncomp.max = 10,
             nbf=0)
}

\arguments{
 \item{X}{a numeric \code{matrix} corresponding to a library of mass spectra. Each row of \code{X} is the intensities of a mass spectrum measured on mass-over-charge values. The columns are assumed to be mass-over-charge values.}

 \item{Y}{a \code{factor} with a length equal to the number of rows in \code{X} and containing the categories of each mass spectrum in \code{X}.}

  \item{MethodSelection}{a \code{character} indicating the method used for variables selection. Six methods are available: (1) \code{"RFERF"} for recursive feature elimination (RFE) coupled with random forests (see \code{rfe} in the \code{caret} R package); (2) \code{"RFEGlmnet"} for RFE with coupled with logistic regression; (3) \code{"VSURF"} for a method using random forests (see \code{VSURF} in the \code{VSURF} R package); (4) \code{"sPLSDA"} for a method based on sparse partial least squares discriminant analysis (see \code{splsda} in the \code{mixOmics}); (5) \code{"mda"} for a method selecting variables from the distribution of the "mean decrease in accuracy" variables importances of a random forest (see \code{importance} function in the \code{randomForest} R package); (6) \code{"cvp"} for a method selecting variables from the distribution of the cross-validated  permutation variables importances of a random forest (see \code{CVPVI} function in the \code{vita} R package). Additional explanations are available in the Details section.}

  \item{MethodValidation}{a \code{character} indicating the resampling method:\code{"cv"} for cross-validation; \code{"repeatedcv"} for repeated cross-validation; and \code{"LOOCV"} for leave-one-out cross-validation. Only used for the \code{"RFERF"}, \code{"RFEGlmnet"} and \code{"sPLSDA"} methods.}

  \item{NumberCV}{a \code{numeric} value indicating the number of K-folds for cross-validation. Only used for the \code{"RFERF"}, \code{"RFEGlmnet"}, \code{"sPLSDA"} and \code{"cvp"} methods.}

   \item{RepeatsCV}{a \code{numeric} value indication the number of repeat(s) for K-folds for cross-validation or repeated cross-validation. Only used for the \code{"RFERF"}, \code{"RFEGlmnet"} and \code{"sPLSDA"} methods.}

   \item{PreProcessing}{a \code{vector} indicating the method(s) used to pre-process the mass spectra in \code{X}: centering (\code{"center"}), scaling (\code{"scale"}), eliminating near zero variance predictors (\code{"nzv"}), or correlated predictors (\code{"corr"}). Only used for the \code{"RFERF"}, \code{"RFEGlmnet"} and \code{"sPLSDA"} methods.}
   
   \item{Metric}{a \code{character} indicating the metric used to select the optimal model for the RFE algorithms. Possible metrics are the \code{"Kappa"} coefficient or the \code{"Accuracy"}. This argument is not used for the \code{"VSURF"}, \code{"cvp"}, \code{"mda"} and the \code{"sPLSDA"} methods of \code{MethodSelection}. See details of the \code{"SelectionVar"} function.}
   
    \item{Sampling}{a \code{character} indicating an optional subsampling method to handle imbalanced datasets: subsampling methods are either \code{"no"} (no subsampling), \code{"up"}, \code{"down"} or \code{"smote"}. \code{"no"} by default.}
    
  \item{Sizes}{a numeric \code{vector} indicating the number of variables to select. Only used for the \code{"RFERF"}, \code{"RFEGlmnet"} and \code{"sPLSDA"} methods. For the \code{"RFERF"} and \code{"RFEGlmnet"} methods, the final number of selected variables is the one giving the highest average \code{"Metric"} (\code{"Accuracy"} or \code{"Kappa"}) on the folds used for cross-validation. It is thus bounded by \code{NumberCV*max(Sizes)}. For the \code{"sPLSDA"} method, \code{Sizes} corresponds to the number of variables to test from the \code{X} dataset when estimating the sparse PLS-DA model (see \code{test.keepX} argument in the \code{mixOmics} R package).}
  
  \item{Ntree}{a \code{numeric} value indicating the number of trees in random forests, only used if \code{MethodSelection} = \code{"VSURF"} or \code{"mda"} or \code{"cvp"}. Note we advise to select a number highly superior to the total number of variables for a robust selection (to not miss some features in the subspaces used to build trees). It is 1000 by default.}
  
  \item{ncores}{a \code{positive integer} only used for the \code{cvp} method. The number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores. If \code{ncores=0}, then the half of CPU cores on the current host are used.}

  \item{ncomp.max}{a \code{positive integer} indicating the maximum number of components that can be included in the sPLS-DA model (10 by default).}

  \item{threshold}{a \code{numeric} value corresponding to a threshold used for the optimal selection of the number of components included in the sPLS-DA model (0.01 by default). When the number of components increases and the balanced classification error rate (BER) does not change anymore, we keep the minimal number where the BER reaches a plateau (i.e. when \code{BER(N)-BER(N+1)<threshold}, we keep \code{N}). If a plateau is not reached, \code{ncomp.max} components are selected.}
  
  \item{nbf}{a \code{numeric} value corresponding to a number of simulated non discriminant features. This is used to improve the robustness of the estimation of the distribution of the variable importances for non discriminant features. Only used for the \code{"mda"} and \code{"cvp"} methods. 0 by default: no additional non discriminant feature is created.}
}

\value{A list composed of:

\item{sel_moz}{a \code{vector} with discriminant mass-over-chage values.}

For the \code{"RFERF"} and \code{"RFEGlmnet"} methods, it also returns the results of the \code{rfe} function of the \code{caret} R package.

For the \code{"VSURF"} method, it also returns the results of the results of the \code{VSURF} function of the \code{VSURF} R package. 

For the \code{"sPLSDA"} method, it also returns the following items:
  \item{Raw_data}{a horizontal bar plot and  containing the contribution of features on each component.}

  \item{selected_variables}{\code{data frame} with uniques features (selected variables to keep and containing the contribution of features in order to class samples).See \code{plotLoadings} in the \code{mixOmics} R package for details.}

For the \code{"mda"} and \code{"cvp"} methods, it also returns the following items:

 \item{nb_to_sel}{a numeric value corresponding to an estimated number of mass-over-chage values where the intensities are significantly different between categories (see details).}
 
 \item{imp_sel}{a vector containing the variable importances for the selected features.}
}
\details{
The selection of variables can be carried out with two different objectives: either to find a minimum number of variables allowing to obtain the highest possible accuracy (or Kappa coefficient), which involves the possible elimination of variables correlated between them (i.e. not bringing any additional predictive power with respect to some other variables); or to find all the variables in the dataset with a potential predictive power ("discriminant" variables).

The \code{VSURF} method attempts to accomplish only the first objective. 
The \code{mda} and \code{cvp} methods attempt to accomplish the second objective, as do the methods available in the \code{SelectionVarStat} function of our \code{MSclassifR} R package.
The \code{RFERF}, \code{RFEGlmnet} and \code{sPLSDA} methods take as input a number of variables to be selected(\code{Sizes} argument), and can therefore be used with both objectives.

Within the framework of the second objective, either the \code{mda} or \code{cvp} methods can be used to estimate a number of discriminant variables from the importances of variables. The \code{SelectionVarStat} function can also be used to estimate this number from distributions of p-values. Of note, be sure that the \code{Ntree} argument is high enough to get a robust estimation with the \code{mda} or \code{cvp} methods.

The \code{"RFEGlmnet"} and \code{"RFERF"} methods are based on recursive feature elimination and can either optimize the kappa coefficient or the accuracy as metrics when selecting variables.

The \code{"sPLSDA"} method selects variables from the ones kept in latent components of the sparse PLS-DA model using an automatic choice of the number of components (when the balanced classification error rate (BER) reaches a plateau - see argument \code{threshold}).

The \code{"mda"} and \code{"cvp"} methods use the distribution of variable importances to estimate the number of discriminant features (mass-over-charge values). Briefly, the distribution of variable importances for useless (not discriminant) features is firstly estimated from negative importance variables by the method proposed in section 2.6 of Janitza et al.(2018). Next, the following mixture model is assumed:
\eqn{F(x)=\pi\times F_u(x)+(1-\pi)\times F_d(x)}
where \eqn{F} is the empirical cumulative distribution of variable importances of all the features, \eqn{F_u} the one of the useless features, \eqn{F_d} the one of the discriminative features, and \eqn{\pi} is the proportion of useless features in the dataset. 
From the estimated distribution of useless features, we can estimate quantile values \eqn{x_q} and compute \eqn{\epsilon_q=min(F(x_q)/q;1)} for each quantile \eqn{q}. The minimum of the \eqn{\epsilon_q} corresponds to the estimated proportion of useless features in the dataset, what allows estimating the number of discriminant features by \eqn{N_d=floor(N\times (1 - \pi))} where N is the total number of features. Next, the \eqn{N_d} features with the highest variable importances are selected.

The \code{"VSURF"} and \code{"sPLSDA"} methods use the minimum mean out-of-bag (OOB) and balanced classification error rate (BER) metrics respectively. 

For \code{Sampling} methods available for unbalanced data: \code{"up"} corresponds to the up-sampling method which consists of random sampling (with replacement) so that the minority class is the same size as the majority class; \code{"down"} corresponds to the down-sampling method randomly which consists of random sampling (without replacement) of the majority class so that their class frequencies match the minority class; \code{"smote"} corresponds to the Synthetic Minority Over sampling Technique (SMOTE) specific algorithm for data augmentation which consist of creates new data from minority class using the K Nearest Neighbor algorithm.

See \code{rfe} in the \code{caret} R package, \code{VSURF} in the \code{VSURF} R package, \code{splsda} in the \code{mixOmics} R package, \code{importance} function in the \code{randomForest} R package, and \code{CVPVI} function in the \code{vita} R package for more details.
}

\references{

Kuhn, Max. (2012). The caret Package. Journal of Statistical Software. 28.

Genuer, Robin, Jean-Michel Poggi and Christine Tuleau-Malot. VSURF : An R Package for Variable Selection Using Random Forests. R J. 7 (2015): 19.

Friedman J, Hastie T, Tibshirani R (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22.

Kim-Anh Le Cao, Florian Rohart, Ignacio Gonzalez, Sebastien Dejean with key contributors Benoit Gautier, Francois, Bartolo, contributions from Pierre Monget, Jeff Coquery, FangZou Yao and Benoit Liquet. (2016). mixOmics: Omics. Data Integration Project. R package version 6.1.1. \url{https://CRAN.R-project.org/package=mixOmics}

Nitesh V. Chawla, Kevin W. Bowyer, Lawrence O. Hall, and W. Philip Kegelmeyer. 2002. SMOTE: synthetic minority over-sampling technique. J. Artif. Int. Res. 16, 1 (January 2002), 321–357.

Branco P, Ribeiro R, Torgo L (2016). “UBL: an R Package for Utility-Based Learning.” CoRR, abs/1604.08079. 

Janitza, S., Celik, E., Boulesteix, A. L. (2018). A computationally fast variable importance test for random forests for high-dimensional data. Advances in Data Analysis and Classification, 12, 885-915.
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
b <- SelectionVar(X, Y, MethodSelection = c("VSURF"))
summary(b$result)

###############################################################################
## 4. Perform variables selection using SelectionVar with "mda" or "cvp"
# option 1: Using mean decrease in accuracy  
# with no sampling method
c <- SelectionVar(X,Y,MethodSelection="mda",Ntree=10*ncol(X)) 

# Estimation of the number of peaks to discriminate species
c$nb_to_sel

# Discriminant mass-over-charge values 
c$sel_moz

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
            Peaks2=c$sel_moz,col_spec="blue",col_peak="black")

# option 2: Using cross-validated permutation variable importance measures (more "time-consuming")       
# with no sampling method
d <- SelectionVar(X,Y,MethodSelection="cvp",NumberCV=2,ncores=2,Ntree=1000)

# Estimation of the number of peaks to discriminate species
d$nb_to_sel

# Discriminant mass-over-charge values 
d$sel_moz

# Plotting peaks on the first pre-processed mass spectrum and highlighting the 
# discriminant mass-over-charge values with red lines
PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],
            Peaks2=d$sel_moz,col_spec="blue",col_peak="black")

# Mass-over charge values found with both methods ("mda" and "cvp")
intersect(c$sel_moz,d$sel_moz)


}

}

\keyword{Feature selection}
