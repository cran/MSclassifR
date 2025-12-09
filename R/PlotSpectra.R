#' Plot spectral data with optional peak markers
#'
#' Create a ggplot of a mass spectrum with optional points for detected peaks
#' and optional vertical lines (and points) highlighting user-specified m/z values.
#'
#' @param SpectralData An object containing spectrum data with numeric slots
#'   `@mass` (m/z) and `@intensity` (signal). Typically a MALDIquant
#'   MassSpectrum-like S4 object.
#' @param absx Either the string "ALL" to plot the full m/z range, or a numeric
#'   length-2 vector c(min, max) specifying the m/z window to display.
#' @param Peaks Optional MassPeaks object (e.g., from MALDIquant) providing
#'   detected peak positions (`@mass`) and intensities (`@intensity`) to plot as points.
#' @param Peaks2 Optional numeric or character vector of m/z values to highlight.
#'   If character/factor, values are coerced to numeric (non-numeric chars removed).
#'   Vertical dashed lines are drawn at these m/z. If `Peaks` is also supplied,
#'   points are plotted for peaks whose m/z match `Peaks2` (within `tol`).
#' @param col_spec Colour for the spectrum line. Default: 1 (black).
#' @param col_peak Colour for points corresponding to all peaks in `Peaks`.
#'   Default: 2 (red).
#' @param shape_peak Point shape for `Peaks` points. Default: 3.
#' @param col_peak2 Colour for points corresponding to the subset of `Peaks`
#'   that match `Peaks2`. Default: 2 (red).
#' @param shape_peak2 Point shape for the `Peaks2`-matched points. Default: 2.
#' @param tol Numeric tolerance (in m/z units) used to match `Peaks@mass` to
#'   `Peaks2`. Set to 0 for exact matching (may miss due to floating-point
#'   precision). Default: 0.
#'
#' @return A ggplot object representing the spectrum and optional annotations.
#' @examples
#' \donttest{
#' if (requireNamespace("MALDIquant", quietly = TRUE)) {
#' # Load mass spectra
#'  data("CitrobacterRKIspectra", package = "MSclassifR")
#' # Plot raw mass spectrum
#'  PlotSpectra(SpectralData = CitrobacterRKIspectra[[1]])
#' # standard pre-processing of mass spectra
#'  spectra <- SignalProcessing(CitrobacterRKIspectra)
#' # Plot pre-processed mass spectrum
#'  PlotSpectra(SpectralData=spectra[[1]])
#' # detection of peaks in pre-processed mass spectra
#'  peaks <- PeakDetection(x = spectra, averageMassSpec=FALSE)
#' # Plot peaks on pre-processed mass spectrum
#'  PlotSpectra(SpectralData=spectra[[1]],Peaks=peaks[[1]],col_spec="blue",col_peak="black")
#' }
#' }
#' @export
PlotSpectra <- function(SpectralData,
                        absx = "ALL",
                        Peaks = NULL,
                        Peaks2 = NULL,
                        col_spec = 1,
                        col_peak = 2,
                        shape_peak = 3,
                        col_peak2 = 2,
                        shape_peak2 = 2,
                        tol = 0) {
  # for NSE notes
  mass <- intensity <- peak <- NULL

  # Validate/derive plotting range
  if (is.character(absx) && length(absx) == 1 && toupper(absx) == "ALL") {
    x_range <- range(SpectralData@mass, na.rm = TRUE)
  } else if (is.numeric(absx) && length(absx) == 2) {
    x_range <- range(absx)
  } else {
    stop("absx must be 'ALL' or a numeric vector of length 2.")
  }

  # Base spectrum DF (filtered by range)
  DF <- data.frame(mass = SpectralData@mass, intensity = SpectralData@intensity)
  DF <- DF[DF$mass >= x_range[1] & DF$mass <= x_range[2], , drop = FALSE]

  # Peaks points (all detected peaks)
  DFpeaks <- NULL
  if (!is.null(Peaks)) {
    if (!inherits(Peaks, "MassPeaks")) stop("'Peaks' must be a 'MassPeaks' object.")
    DFpeaks <- data.frame(mass = Peaks@mass, peak = Peaks@intensity)
    DFpeaks <- DFpeaks[DFpeaks$mass >= x_range[1] & DFpeaks$mass <= x_range[2], , drop = FALSE]
  }

  # Peaks2: specific m/z to highlight (accept numeric/character/factor)
  DFpeaks2 <- NULL
  lines_mz <- NULL
  if (!is.null(Peaks2)) {
    # Coerce Peaks2 to numeric if needed
    if (is.factor(Peaks2)) Peaks2 <- as.character(Peaks2)
    if (!is.numeric(Peaks2)) {
      # first try direct as.numeric, then strip non-numeric chars
      p2_try <- suppressWarnings(as.numeric(Peaks2))
      if (any(is.na(p2_try))) {
        p2_try2 <- suppressWarnings(as.numeric(gsub("[^0-9\\.]+", "", as.character(Peaks2))))
        Peaks2 <- p2_try2
      } else {
        Peaks2 <- p2_try
      }
    }
    # After coercion, keep finite values only
    Peaks2 <- Peaks2[is.finite(Peaks2)]
    # limit vertical lines to range
    lines_mz <- Peaks2[Peaks2 >= x_range[1] & Peaks2 <= x_range[2]]

    # If Peaks provided, derive points for Peaks that match Peaks2 (within tol)
    if (!is.null(Peaks) && length(lines_mz)) {
      if (tol > 0) {
        keep <- rowSums(abs(outer(Peaks@mass, lines_mz, "-")) <= tol) > 0
        if (any(keep)) {
          DFpeaks2 <- data.frame(mass = Peaks@mass[keep], peak = Peaks@intensity[keep])
          DFpeaks2 <- DFpeaks2[DFpeaks2$mass >= x_range[1] & DFpeaks2$mass <= x_range[2], , drop = FALSE]
        }
      } else {
        idx <- Peaks@mass %in% lines_mz
        if (any(idx)) {
          DFpeaks2 <- data.frame(mass = Peaks@mass[idx], peak = Peaks@intensity[idx])
        }
      }
    }
  }

  # Plot
  p <- ggplot2::ggplot(DF, ggplot2::aes(x = mass, y = intensity)) +
    ggplot2::geom_line(colour = col_spec) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "mass over charge (m/z)")

  if (!is.null(DFpeaks) && nrow(DFpeaks)) {
    p <- p + ggplot2::geom_point(
      data = DFpeaks,
      mapping = ggplot2::aes(x = mass, y = peak),
      colour = col_peak,
      shape = shape_peak
    )
  }

  if (!is.null(DFpeaks2) && nrow(DFpeaks2)) {
    p <- p + ggplot2::geom_point(
      data = DFpeaks2,
      mapping = ggplot2::aes(x = mass, y = peak),
      colour = col_peak2,
      shape = shape_peak2
    )
  }

  if (!is.null(lines_mz) && length(lines_mz)) {
    # No aes() here -> avoids NOTE: no visible binding for 'xint'
    p <- p + ggplot2::geom_vline(
      xintercept = lines_mz,
      linetype = 2, colour = 2, linewidth = 0.5
    )
  }

  return(p)
}
