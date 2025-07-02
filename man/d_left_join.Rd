\name{d_left_join}

\alias{d_left_join}

\title{Function joining two tables based not on exact matches}

\description{
This function joins two tables based on a distance metric of one or more columns. It gives the same results that the \code{distance_left_join} function of the \code{fuzzyjoin} R package, but its execution time is faster. It is used to match mass spectra to short-listed mass-to-charge values in \code{PredictLogReg} and \code{PredictFastClass} functions.
}

\usage{

d_left_join (x, y, by = NULL, method = "euclidean", max_dist = 1,
             distance_col = NULL)

}

\arguments{
  \item{x}{A \code{data.frame} object.}

  \item{y}{A \code{data.frame} object.}

  \item{by}{Columns by which to join the two tables. \code{NULL} by default (tables are joined by matching the common column names that appear in both datasets x and y).}

  \item{method}{Method to use for computing distance, either "euclidean" (default) or "manhattan".}

  \item{max_dist}{A \code{numeric} value indicating the maximum distance to use for joining. 1 by default.}

  \item{distance_col}{A \code{character} that specifies the name of a new column to be added to the output, which will contain the calculated distances   between both tables. If NULL, no column is added (default).}
}

\value{Returns a data frame that contains the results of joining the \code{x} dataset onto the \code{y} datasets, where all rows from \code{x} are preserved and matching data from \code{y} is added according to the specified criteria.}






