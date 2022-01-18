##  Plot spectral data ##


PlotSpectra <- function (SpectralData, absx="ALL", Peaks=NULL, Peaks2=NULL, col_spec=1, col_peak=2,
                         shape_peak=3, col_peak2=2, shape_peak2=2){
  peak=NULL
  ## Create a data frame
  if (length(absx)==1){
  if (absx=="ALL"){

      DF <- data.frame(mass = SpectralData@mass, intensity = SpectralData@intensity)

      if (!is.null(Peaks)){
        if (class(Peaks)=="MassPeaks"){
      DFpeaks=data.frame(mass=Peaks@mass,peak=Peaks@intensity)
        }}

      if (!is.null(Peaks2)){
        if (is.vector(Peaks2)){
      DFpeaks2=data.frame(mass=Peaks@mass[which(Peaks@mass%in%Peaks2)],peak=Peaks@intensity[which(Peaks@mass%in%Peaks2)])
        }}

  }}else{

      DF <- data.frame(mass = SpectralData@mass , intensity = SpectralData@intensity )
      sel=which((DF$mass>min(absx))&(DF$mass<max(absx)))
      DF = DF [sel,]
      #
      if (!is.null(Peaks)){
        if (class(Peaks)=="MassPeaks"){
          DFpeaks=data.frame(mass=Peaks@mass,peak=Peaks@intensity)
          sel=which((DFpeaks$mass>min(absx))&(DFpeaks$mass<max(absx)))
          DFpeaks = DFpeaks [sel,]
          #
        }}

      if (!is.null(Peaks2)){
        if (is.vector(Peaks2)){
          DFpeaks2=data.frame(mass=Peaks@mass[which(Peaks@mass%in%Peaks2)],peak=Peaks@intensity[which(Peaks@mass%in%Peaks2)])
          sel=which((DFpeaks2$mass>min(absx))&(DFpeaks2$mass<max(absx)))
          DFpeaks2 = DFpeaks2 [sel,]
          #
        }}
      #
      Peaks2 = Peaks2 [which((Peaks2>min(absx))&(Peaks2<max(absx)))]
  }

  ## Plot with ggplot2
  mass <- NULL
  intensity <- NULL
  spectra_plot <- ggplot2::ggplot(DF, ggplot2::aes(x = mass, y = intensity))
  spectra_plot <- spectra_plot + ggplot2::geom_line(col=col_spec) +  ggplot2::theme_bw()
  spectra_plot <- spectra_plot + xlab("mass over charge (m/z)")

  if (!is.null(Peaks)){
     if (class(Peaks)=="MassPeaks"){
        spectra_plot <- spectra_plot + ggplot2::geom_point(data=DFpeaks,
                                            mapping=ggplot2::aes(x = mass, y = peak),col=col_peak,shape=shape_peak)
     }}

  if (!is.null(Peaks2)){
      if (is.vector(Peaks2)){
          spectra_plot <- spectra_plot + ggplot2::geom_point(data=DFpeaks2,
                                                         mapping=ggplot2::aes(x = mass, y = peak),col=col_peak2,shape=shape_peak2)
          for (i in 1:length(Peaks2)){
              spectra_plot <- spectra_plot + ggplot2::geom_vline(xintercept=as.numeric(Peaks2)[i],
                                                                   linetype=2, color=2, size=0.5)
          }
      }}

  return(spectra_plot)
}




