\name{calculate_distance}

\alias{calculate_distance}

\title{Function calculating the distance between two vectors.}

\description{
This function calculates the distance between two vectors using the specified distance metric. Distance metrics available are "p-norm", "Chebyshev", "Canberra", "Overlap", "HEOM" (Heterogeneous Euclidean-Overlap Metric), "HVDM" (Heterogeneous Value Difference Metric). Used in the \code{fast_find_neighbors} function of our package.
}

\usage{

calculate_distance(x, y, nominal_indices, p_code)

}

\arguments{
  \item{x}{Vector of numeric and/or categorical values.}

  \item{y}{Vector of numeric and/or categorical values.}

  \item{nominal_indices}{Vector indicating which positions in \code{x} and \code{y} contain categorical variables.This distinction is needed because:
   \itemize{
     \item For categorical variables, distances are often based on whether values match or not.
     \item Hybrid distance metrics like HEOM and HVDM require knowing which variables are nominal to apply appropriate distance calculations.
  }}

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

\value{A numeric value representing the distance between the two vectors.}

\details{Different distance metrics handle nominal variables differently:
 \itemize{
   \item For pure numeric metrics (p >= 1, p = 0, p = -1), nominal features are ignored
   \item For the Overlap metric (p = -2), only nominal features are considered
   \item For HEOM (p = -3), numeric features use normalized Euclidean distance while
     nominal features use overlap distance (1 if different, 0 if same)
   \item For HVDM (p = -4), a specialized metric combines normalized differences for
     numeric features and value difference metric for nominal features
}}





