#' Internal helpers and C++ exports
#'
#' These functions are internal building blocks used by SignalProcessingUltra:
#'
#' - alignSpectra_landmark_cpp: wrapper that applies landmark matching/warping per spectrum.
#' - cpp_match_anchors: Rcpp function for greedy 1:1 matching between spectrum peaks and reference anchors.
#' - cpp_warp_mz_piecewise: Rcpp function for piecewise-linear warping of a mass axis.
#'
#' They are not part of the user-facing API.
#' @keywords internal
#' @name msclassifr-internal
#' @aliases alignSpectra_landmark_cpp cpp_match_anchors cpp_warp_mz_piecewise
NULL
