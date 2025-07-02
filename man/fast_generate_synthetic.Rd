\name{fast_generate_synthetic}

\alias{fast_generate_synthetic}

\title{Function generating synthetic examples using SMOTE}

\description{
This function generates synthetic examples using the SMOTE algorithm.
For each example in the minority class, this function generates synthetic examples by interpolating between the example and its k nearest neighbors. This function is used in the \code{smote_classif} function.}

\usage{

fast_generate_synthetic(dat, k, n, p_code)

}

\arguments{
  \item{dat}{A \code{data.frame} object.}

  \item{k}{Number of nearest neighbors to consider.}

  \item{n}{Number of synthetic examples to generate.}

  \item{p_code}{Numeric code representing the distance metric to use:
    \itemize{
     \item p >= 1: p-norm
     \item p = 0: Chebyshev
     \item p = -1: Canberra
     \item p = -2: Overlap (nominal attributes only)
     \item p = -3: HEOM (Heterogeneous Euclidean-Overlap Metric)
     \item p = -4: HVDM (Heterogeneous Value Difference Metric)
   }}
}

\value{A data frame containing the generated synthetic examples.}

\seealso{\code{\link{smote_classif}}}




