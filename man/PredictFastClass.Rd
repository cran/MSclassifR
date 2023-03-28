\name{PredictFastClass}

\alias{PredictFastClass}

\title{Prediction of the category to which a mass spectrum belongs using linear regressions of mass spectra.}

\description{For each mass peak in a list of mass peaks, a linear regression is performed between the mass spectrum and mass spectra corresponding to a category. This is performed for each category and associated to an Akaike Information Criterium.  Next, the AIC are used to determine the belonging of a mass spectrum to a category. It also provides a probability that the mass spectrum does not belong to any of the input categories.  
}


\usage{
PredictFastClass(peaks, 
                 mod_peaks,
                 Y_mod_peaks,
                 moz="ALL",
                 tolerance = 6,
                 toleranceStep = 2,
                 normalizeFun = TRUE,
                 noMatch = 0)  
}

\arguments{
  \item{peaks}{a list of \code{MassPeaks} objects (see \code{MALDIquant} R package).}
  \item{mod_peaks}{an intensity matrix corresponding to mass spectra for which the category is known. Each column is a mass-over-charge value, each row corresponds to a mass spectrum.}
  \item{Y_mod_peaks}{a \code{factor} with a length equal to the number of mass spectra in \code{mod_peaks} and containing the categories of each mass spectrum in \code{mod_peaks}.}
  \item{moz}{a \code{vector} with the set of shortlisted mass-over-charge values that corresponds to mass-over-charge values in the columns of \code{mod_peaks}. By default, all the mass-over-charge values in \code{mod_peaks} are used.}
  \item{tolerance}{a \code{numeric} value of accepted tolerance to match peaks to the set of shortlisted mass-over-charge values. It is fixed to 6 Da by default.}
  \item{toleranceStep}{a \code{numeric} value added to the \code{tolerance} parameter to match peaks to the set of shortlisted mass-over-charge values. It is fixed to 2 Da by default.}
  \item{normalizeFun}{a \code{logical} value, if \code{TRUE} (default) the maximum intensity will be equal to 1, the other intensities will be expressed in ratio to this maximum.}
  \item{noMatch}{a \code{numeric} value used to replace intensity values if there is no match detected between peaks and the set of shortlisted mass-over-charge values \code{moz}. It is fixed to 0 by default.}
}

\value{Returns a \code{dataframe} containing AIC criteria by category for each mass spectrum in \code{peaks}. The AIC criterion should be minimal for the most probable category. The \code{pred_cat} column is the predicted category for each mass spectrum in \code{peaks}. The \code{p_not_in_DB} is the minimal p-value of several Fisher tests testing if all the linear coefficients associated to mass spectra of a category are null. It can be interpreted as a p-value that the mass spectrum is not present in the input database.}


\examples{
\donttest{

library("MSclassifR")
library("MALDIquant")

# load mass spectra and their metadata
data("CitrobacterRKIspectra","CitrobacterRKImetadata", package = "MSclassifR")
# standard pre-processing of mass spectra
spectra <- SignalProcessing(CitrobacterRKIspectra)
# detection of peaks in pre-processed mass spectra
peaks <- peaks <- MSclassifR::PeakDetection(x = spectra, averageMassSpec=FALSE)
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

#Predict species without peak selection using a tolerance of 1 Da
res = PredictFastClass(peaks=peaks[1:5],
                       mod_peaks=X,
                       Y_mod_peaks=Y,
                       tolerance = 1)

#comparing predicted categories (species) and the truth
cbind(res$pred_cat,as.character(Y[1:5]))

# The method can be applied after a peak selection step
a <- SelectionVar(X,
                  Y,
                  MethodSelection = c("RFERF"),
                  MethodValidation = c("cv"),
                  PreProcessing = c("center","scale","nzv","corr"),
                  NumberCV = 2,
                  Metric = "Kappa",
                  Sizes = c(20:40),
                  Sampling = "up")

#Predict species from selected peaks using a tolerance of 1 Da
res = PredictFastClass(peaks=peaks[1:5],
                       moz = a$sel_moz,
                       mod_peaks=X,
                       Y_mod_peaks=Y, tolerance = 1)

#comparing predicted categories (species) and the truth
cbind(res$pred_cat,as.character(Y[1:5]))

}
}

\keyword{Prediction}

