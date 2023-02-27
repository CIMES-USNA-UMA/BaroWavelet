PlotHRV <-
  function(hrv,
           time_flags = NULL,
           col = "brown",
           tem  = FALSE,
           plotHF = TRUE,
           plotLF = TRUE) {
    Time <- hrv$Time
    HF <- hrv$HF
    LF <- hrv$LF
    LFHF <- hrv$LFHF
    if (tem) {
      im <- tempfile(fileext = ".png")
      png(
        filename = im,
        width = 6,
        height = 6,
        units = "in",
        res = 400
      )
    }
    plottingData <-
      data.frame(
        Time = Time,
        LF = LF,
        HF = HF,
        LFHF = LFHF
      )
    if (plotHF & !plotLF) {
      im <-
        ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = HF)) +
        ggplot2::geom_line()
      return(im)
    } else if (plotLF & !plotHF) {
      im <-
        ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = LF)) +
        ggplot2::geom_line()
      return(im)
    } else if (plotLF & plotHF) {
      im <-
        ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = LFHF)) +
        ggplot2::geom_line()
      return(im)
    }
    
    if (tem) {
      dev.off()
      return(im)
    }
    
  }
