\name{fast_find_neighbors}

\alias{fast_find_neighbors}

\title{Function finding k Nearest Neighbors for each row of a matrix}

\description{
This function finds the k nearest neighbors for each row in a matrix using the specified distance metric. Distance metrics available are "p-norm", "Chebyshev", "Canberra", "Overlap", "HEOM" (Heterogeneous Euclidean-Overlap Metric), "HVDM" (Heterogeneous Value Difference Metric). See the \code{calculate_distance} function of our package for more details on the distances.
}

\usage{

fast_find_neighbors(data, nominal_indices, p_code, k)

}

\arguments{
  \item{data}{Matrix where the k nearest neighbors for each row are searched.}

  \item{nominal_indices}{Vector of column indices indicating which features are
   categorical (nominal) variables. This is crucial for proper distance calculation
   as nominal and numeric features require different handling. For example, if columns
   2 and 5 contain categorical variables, nominal_indices should be \code{c(2, 5)}. See \code{calculate_distance} function.}

  \item{p_code}{Numeric code representing the distance metric to use:
    \itemize{
     \item p >= 1: p-norm
     \item p = 0: Chebyshev
     \item p = -1: Canberra
     \item p = -2: Overlap (nominal attributes only)
     \item p = -3: HEOM (Heterogeneous Euclidean-Overlap Metric)
     \item p = -4: HVDM (Heterogeneous Value Difference Metric)
   }}

  \item{k}{Number of nearest neighbors to find.}

}

\value{A matrix where each row contains the indices of the k nearest neighbors for the corresponding example.}

\seealso{\code{\link{calculate_distance}}}




