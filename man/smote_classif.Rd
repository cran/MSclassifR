\name{smote_classif}

\alias{smote_classif}

\title{SMOTE for Classification Problems}

\description{
This function performs Synthetic Minority Over-sampling Technique (SMOTE) to address class imbalance in classification problems. This implementation supports various distance metrics, different balancing strategies, and handles mixed data types (numeric and categorical). Its execution time is faster than the \code{SmoteClassif} of the \code{UBL} R package.}

\usage{

smote_classif (form, dat, C.perc = "balance", k = 5, repl = FALSE,
                          dist = "Euclidean", p = 2)

}

\arguments{
  \item{form}{A model formula identifying the target variable (e.g., \code{Class ~ .}).}

  \item{dat}{A data frame containing the imbalanced dataset.}

  \item{C.perc}{Either \code{"balance"}, \code{"extreme"}, or a named list containing over/under-sampling
percentages for each class. Values < 1 indicate under-sampling, values > 1 indicate over-sampling, and values = 1 indicate no change. "balance" equalizes all classes, \code{"extreme"} performs more aggressive balancing.}

  \item{k}{Integer specifying the number of nearest neighbors to use when generating synthetic examples (default: 5).}

   \item{repl}{Logical, whether to allow sampling with replacement when under-sampling (default: FALSE).}

   \item{dist}{Distance metric to use for nearest neighbor calculations. Supported metrics: \code{"Euclidean"} (default), \code{"Manhattan"}, \code{"Chebyshev"}, \code{"Canberra"}, \code{"Overlap"}, \code{"HEOM"}, \code{"HVDM"}, or \code{"p-norm"}. See \code{calculate_distance} function for details.}

   \item{p}{Parameter used when \code{dist = "p-norm"} (default: 2).}
}

\value{A data frame with the same structure as the input, but with rebalanced classes according to the specified strategy.}

\examples{

# Load the iris dataset
data(iris)

# Create an imbalanced dataset by taking a subset
imbal_iris <- iris[c(1:40, 51:100, 101:110), ]
table(imbal_iris$Species)  # Show class distribution

# Balance classes using the default "balance" strategy
balanced_iris <- smote_classif(Species ~ ., imbal_iris)
table(balanced_iris$Species)  # Show balanced distribution

# Custom over/under-sampling
custom_iris <- smote_classif(Species ~ ., imbal_iris,
                            C.perc = list(setosa = 0.8,
                                          versicolor = 1,
                                          virginica = 3))
table(custom_iris$Species)

}

\details{If you use this package in your research, please cite the associated publication (\doi{10.1016/j.eswa.2025.128796}).}

\references{Chawla, N. V., Bowyer, K. W., Hall, L. O., & Kegelmeyer, W. P. (2002).
SMOTE: Synthetic Minority Over-sampling Technique. Journal of Artificial
Intelligence Research, 16, 321-357.\doi{10.1613/jair.953}.

Alexandre Godmer, Yahia Benzerara, Emmanuelle Varon, Nicolas Veziris, Karen Druart, Renaud Mozet, Mariette Matondo, Alexandra Aubry, Quentin Giai Gianetto, MSclassifR: An R package for supervised classification of mass spectra with machine learning methods, Expert Systems with Applications, Volume 294, 2025, 128796, ISSN 0957-4174, \doi{10.1016/j.eswa.2025.128796}.
}

\seealso{\code{\link{fast_generate_synthetic}} for the implementation of
synthetic example generation}




