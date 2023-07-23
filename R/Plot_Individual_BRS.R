#' Plot individual BRS
#'
#' Plots an already computed BRS from a specific subject
#' @param fun BRS obtained either by \link[BaroWavelet]{AlphaIndexDWT} or
#' \link[BaroWavelet]{TransferFunCWT}
#' @param avg Boolean. Plot the scale-averaged CWT BRS? Default is FALSE
#' @param time_col Color used to highlight a specific time interval. Default is brown
#' @param HFcolor Color to be used to highlight the HF band for the average CWT plot
#' @param LFcolor Color to be used to highlight the LF band for the average CWT plot
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL
#' @param thr Coherence threshold to be used for the plot. Default is 0.5
#' @param use.thr Boolean, should a coherence threshold be used in the analyses? Default is TRUE.
#' @param nfreqs Number of frequencies to be used. Default is 7
#' @param tem Boolean, creates a temporal file for the plot. Default is FALSE
#' @param newPlot Boolean, generates a new plot without overwriting a previous plot. Default is TRUE
#' @param title Part of the title for the plot. Used to indicate the data source on the title. Default is "data"
#' @param plotHF Boolean, plot results form the HF band. Default is TRUE
#' @param plotLF Boolean, plot results from the LF band. Default is TRUE
#' @param ylim Maximum y axis limit. Default is NULL
#' @param use.ggplot Boolean, use methods from \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2} package to plot the results. Default is FALSE
#' @param HF Maximum limit of the HF band, shown at the plot title. Default is 0.4 Hz
#' @param LF Maximum limit of the LF band, shown at the plot title. Default is 0.15 Hz
#' @param VLF Maximum limit of the VLF band, shown at the plot title. Default is 0.04 Hz
#'
#'
#' @return None
#'
#' @author Alvaro Chao-Ecija
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 labs
#' @importFrom  biwavelet plot.biwavelet
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' AlphaIndex <- AlphaIndexDWT(Data, wv = "d8", error = 0.0005)
#'
#' PlotBRS(AlphaIndex, plotHF = TRUE, plotLF = TRUE, newPlot = FALSE)
PlotBRS <-
  function(fun,
           avg = FALSE,
           time_col = "brown",
           HFcolor = "yellow",
           LFcolor = "green",
           time_flags = NULL,
           thr = 0.5,
           use.thr = TRUE,
           nfreqs = 7,
           tem = FALSE,
           newPlot = TRUE,
           title = "data",
           plotHF = TRUE,
           plotLF = TRUE,
           ylim = NULL,
           use.ggplot = FALSE,
           HF = 0.4,
           LF = 0.15,
           VLF = 0.04) {
    #if(dev.cur() > 1)
    if (newPlot & !tem) {
      dev.new(title = paste("BRS from", title))
    }
    if (fun$type == "brs_dwt") {
      if (!tem & !use.ggplot) {
        PlotDwtBRS(
          fun,
          time_flags = time_flags,
          col = time_col,
          tem = tem,
          plotHF = plotHF,
          plotLF = plotLF,
          ylim = ylim,
          use.ggplot = FALSE,
          fHF = HF,
          fLF = LF,
          fVLF = VLF
        )
      } else {
        im <-
          PlotDwtBRS(
            fun,
            time_flags = time_flags,
            col = time_col,
            tem = tem,
            plotHF = plotHF,
            plotLF = plotLF,
            ylim = ylim,
            use.ggplot = TRUE,
            fHF = HF,
            fLF = LF,
            fVLF = VLF
          )
        return(im)
      }
    } else if (fun$type == "brs_cwt") {
      if (avg) {
        im <- PlotAvgCwtBRS(
          fun,
          thr,
          use.thr,
          time_flags = time_flags,
          len = nfreqs,
          HFcolor = HFcolor,
          LFcolor = LFcolor,
          Tcolor = time_col,
          tem = tem
        )
      } else {
        im <- PlotCwtBRS(fun,
                         thr,
                         use.thr,
                         time_flags = time_flags,
                         tem = tem,
                         Max = ylim)
      }
    }
    if (tem & !newPlot & use.ggplot) {
      return(im)
    }
  }



PlotDwtBRS <-
  function(fun,
           time_flags = NULL,
           col = "brown",
           tem  = FALSE,
           plotHF = TRUE,
           plotLF = TRUE,
           use.xlim = FALSE,
           ylim = NULL,
           use.ggplot = TRUE,
           fHF = 0.4,
           fLF = 0.15,
           fVLF = 0.04) {
    time <- fun$Time
    HF <- fun$HF
    LF <- fun$LF
    if (use.ggplot) {
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
      plottingData <- data.frame(Time = time, LF = LF, HF = HF)
      if (plotHF & plotLF)
        plotHF <- FALSE # Future implementation
      if (plotHF) {
        im <-
          ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = HF)) +
          ggplot2::geom_line()
        if (use.xlim)
          im <-
            im + ggplot2::xlim(time_flags[1] * 60, time_flags[2] * 60)
        if (!is.null(ylim))
          im <- im + ggplot2::ylim(0, ylim)
        #return(im)
      }
      if (plotLF) {
        im <-
          ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = LF)) +
          ggplot2::geom_line()
        if (use.xlim)
          im <-
            im + ggplot2::xlim(time_flags[1] * 60, time_flags[2] * 60)
        if (!is.null(ylim))
          im <- im + ggplot2::ylim(0, ylim)
        #return(im)
      }
      if (tem) {
        dev.off()
        #return(im)
      }
      return(im)
      
    } else {
      if (plotHF & plotLF)
        par(mfrow = c(1, 2))
      if (plotHF) {
        plot(
          time,
          HF,
          type = "l",
          xlab = "Time (s)",
          ylab = "BRS (ms/mmHg)",
          main = paste("HF band (", fHF, " - ", fLF, " Hz)", sep = ""),
          ylim = if (!is.null(ylim))
            c(0, ylim)
        )
        if (is.list(time_flags) &&
            !is.null(col) && (length(time_flags) >= NROW(col))) {
          if (length(time_flags) > NROW(col))
            col <-
              c(col, rep ("black", length(time_flags) - NROW(col)))
          for (ti in 1:length(time_flags)) {
            time_lims <- time_flags[[ti]] * 60
            ti_col <- col[ti]
            limit1 <-
              match(min(abs(time_lims[1] - time)), abs(time_lims[1] - time))
            limit2 <-
              match(min(abs(time_lims[2] - time)), abs(time_lims[2] -  time))
            select_time <- limit1:limit2
            band <- HF
            band[-select_time] <- NA
            lines(time, band, col = ti_col)
          }
          
        } else if (is.numeric(time_flags) &&
                   (NROW(time_flags) == 2) &&
                   !is.null(col) && (NROW(col) == 1)) {
          time_flags <- time_flags * 60
          limit1 <-
            match(min(abs(time_flags[1] - time)), abs(time_flags[1] - time))
          limit2 <-
            match(min(abs(time_flags[2] - time)), abs(time_flags[2] -  time))
          select_time <- limit1:limit2
          band <- HF
          band[-select_time] <- NA
          lines(time, band, col = col)
          time_flags <- time_flags / 60
        }
      }
      if (plotLF) {
        plot(
          time,
          LF,
          type = "l",
          xlab = "Time (s)",
          ylab = "BRS (ms/mmHg)",
          main = paste("LF band (", fLF, " - ", fVLF, " Hz)", sep = ""),
          ylim = if (!is.null(ylim))
            c(0, ylim)
        )
        if (is.list(time_flags) &&
            !is.null(col) && (length(time_flags) >= NROW(col))) {
          if (length(time_flags) > NROW(col))
            col <-
              c(col, rep ("black", length(time_flags) - NROW(col)))
          for (ti in 1:length(time_flags)) {
            time_lims <- time_flags[[ti]] * 60
            ti_col <- col[ti]
            limit1 <-
              match(min(abs(time_lims[1] - time)), abs(time_lims[1] - time))
            limit2 <-
              match(min(abs(time_lims[2] - time)), abs(time_lims[2] -  time))
            select_time <- limit1:limit2
            band <- LF
            band[-select_time] <- NA
            lines(time, band, col = ti_col)
          }
          
        } else if (is.numeric(time_flags) &&
                   (NROW(time_flags) == 2) &&
                   !is.null(col) && (NROW(col) == 1)) {
          time_flags = time_flags * 60
          limit1 <-
            match(min(abs(time_flags[1] - time)), abs(time_flags[1] - time))
          limit2 <-
            match(min(abs(time_flags[2] - time)), abs(time_flags[2] -  time))
          select_time <- limit1:limit2
          band <- LF
          band[-select_time] <- NA
          lines(time, band, col = col)
        }
      }
    }
  }

PlotCwtBRS <-
  function(fun,
           thr = 0.5,
           use.thr = TRUE,
           time_flags = NULL,
           tem  = FALSE,
           Max = NULL,
           use.xlim = FALSE,
           show.coi = TRUE) {
    HF <- fun$HF
    LF <- fun$LF
    VLF <- fun$VLF
    isAlpha <- fun$Alpha
    if (use.thr) {
      mask <- plot.contour <- TRUE
    } else {
      mask <- plot.contour <- FALSE
      thr <- 0
    }
    if (is.null(Max))
      Max <- max(abs(fun$TransferFun))
    fun <- GetBiwaveletObject(fun, thr = thr)
    if (mask) {
      fun$wave[which(fun$rsq < thr)] <- -0.1
    }
    if (plot.contour) {
      col <- c(
        "white",
               "#00007F",
               "blue",
               "#007FFF",
               "cyan",
               "#7FFF7F",
               "yellow",
               "#FF7F00",
               "red",
               "#7F0000"
      )
    } else {
      col <- NULL
    }
    if (tem) {
      im <- tempfile(fileext = ".png")
      png(filename = im,
          width = 1500,
          height = 400)
    }
    xlim <- NULL
    if (use.xlim)
      xlim = time_flags * 60
    par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
    biwavelet::plot.biwavelet(
      fun,
      plot.sig = FALSE,
      plot.phase = TRUE,
      type = "wavelet",
      plot.cb = TRUE,
      arrow.cutoff = thr,
      plot.coi = show.coi,
      main = ifelse(
        isAlpha,
        "Alpha Index by CWT (ms/mmHg)",
        "Transfer Function by CWT (ms/mmHg)"
      ),
      ylim = c(1 / HF, 1 / VLF),
      xlim = xlim,
      fill.cols = col,
      zlim = c(0, Max)
    )
    if (plot.contour) {
      contour(
        fun$t,
        log2(fun$period),
        t(fun$rsq),
        level = thr,
        col = "black",
        lwd = 5,
        add = TRUE,
        drawlabels = FALSE
      )
    }
    if (!is.null(time_flags)) {
      time_flags = time_flags * 60
      abline(v = time_flags[1],
             lwd = 3,
             col = "black")
      abline(v = time_flags[2],
             lwd = 3,
             col = "black")
    }
    abline(h = log2(1 / LF),
           lwd = 3,
           col = "brown")
    if (tem) {
      dev.off()
      return(im)
    }
  }

PlotAvgCwtBRS <- function(fun,
                          thr = 0.5,
                          use.thr = TRUE,
                          scale = 1,
                          time_flags = NULL,
                          len = 7,
                          HFcolor = "yellow",
                          LFcolor = "green",
                          Tcolor = "brown",
                          tem  = FALSE) {
  HF <- fun$HF
  LF <- fun$LF
  VLF <- fun$VLF
  isAlpha <- fun$Alpha
  fun <- GetBiwaveletObject(fun)
  freqs <- 1 / fun$period
  sel_power <- fun$power
  if (use.thr) {
    fun$power[which(fun$rsq < thr)] <- NA
  }
  min_freqs <- min(freqs)
  max_freqs <- max(freqs)
  if (max_freqs > HF)
    max_freqs <- HF
  if (min_freqs < VLF)
    min_freqs <- VLF
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
  par(mfrow = c(3, 1))
  if (!is.null(time_flags)) {
    time_flags <- time_flags * 60
    select_t <-
      fun$t[(fun$t >= time_flags[1]) & (fun$t <= time_flags[2])]
    freq_results <-
      rowMeans(fun$power[, match(select_t, fun$t)], na.rm = TRUE)
  } else {
    freq_results <- rowMeans(fun$power, na.rm = TRUE)
  }
  time_results.HF <-
    colMeans(fun$power[(freqs <= HF) & (freqs > LF),],
             na.rm = TRUE)
  time_results.LF <-
    colMeans(fun$power[(freqs <= LF) & (freqs > VLF),],
             na.rm = TRUE)
  freq_results[is.na(freq_results)] <- 0
  time_results.HF[is.na(time_results.HF)] <- 0
  time_results.LF[is.na(time_results.LF)] <- 0
  plot(
    log2(1 / fun$period),
    freq_results,
    "l",
    xlab = "Frequency",
    ylab = "BRS (ms/mmHg)",
    main = ifelse(
      isAlpha,
      "Frequency Domain Alpha Index",
      "Frequency Domain Transfer Function"
    ),
    ylim =
      c(0, max(freq_results)),
    xaxt = "n"
  )
  polygon(
    x = c(log2(freqs[(freqs <= HF) & (freqs >= VLF)]),
          rev(log2(freqs[(freqs <= HF) & (freqs >= VLF)]))),
    y = c(freq_results[(freqs <= HF) & (freqs >= VLF)],
          double(NROW(freqs[(freqs <= HF) &
                              (freqs >= VLF)]))),
    col = HFcolor
  )
  polygon(
    x = c(log2(freqs[(freqs <= LF) & (freqs > VLF)]),
          rev(log2(freqs[(freqs <= LF) & (freqs > VLF)]))),
    y = c(freq_results[(freqs <= LF) & (freqs > VLF)],
          double(NROW(freqs[(freqs <= LF) & (freqs > VLF)]))),
    col = LFcolor
  )
  axis(1,
       at = log2(1 / fun$period[seq(1, NROW(fun$period), len = len)]),
       labels = round(1 / fun$period[seq(1, NROW(fun$period), len = len)], 3))
  plot(
    fun$t,
    time_results.HF,
    "l",
    xlab = "Time",
    ylab = "BRS (ms/mmHg)",
    main = ifelse(
      isAlpha,
      "Time Domain Alpha Index (HF band)",
      "Time Domain Transfer Function (HF band)"
    )
  )
  if (!is.null(time_flags)) {
    polygon(
      x = c(select_t, rev(select_t)),
      y = c(time_results.HF[match(select_t, fun$t)],
            double(NROW(select_t))),
      col = Tcolor
    )
  }
  plot(
    fun$t,
    time_results.LF,
    "l",
    xlab = "Time",
    ylab = "BRS (ms/mmHg)",
    main = ifelse(
      isAlpha,
      "Time Domain Alpha Index (LF band)",
      "Time Domain Transfer Function (LF band)"
    )
  )
  if (!is.null(time_flags)) {
    polygon(
      x = c(select_t, rev(select_t)),
      y = c(time_results.LF[match(select_t, fun$t)],
            double(NROW(select_t))),
      col = Tcolor
    )
  }
  if (tem) {
    dev.off()
    return(im)
  }
  
}

#########################################################
# OTHER AUXILIARY FUNCTIONS                             #
#########################################################


# Private function: generates an object of class "biwavelet" with a structure
# compatible with the biwavelet package. This is done so that biwavelet methods
# can be applied to the plots.
GetBiwaveletObject <- function(data, use.thr = TRUE, thr = 0.5) {
  biwave_object <- list()
  biwave_object$type <- "wtc"
  biwave_object$t <- data$Time
  biwave_object$period <- 1 / data$Freqs
  biwave_object$power <- abs(data$TransferFun)
  biwave_object$wave <- abs(data$TransferFun)
  #biwave_object$phase <- atan2(Im(data$TransferFun),
  #   Re(data$TransferFun))
  biwave_object$phase <- data$Phase
  biwave_object$rsq <- data$Coherence
  biwave_object$coi <- data$Cone
  if (use.thr) {
    biwave_object$signif <- 5 * (biwave_object$rsq >= thr)
  }
  class(biwave_object) <- "biwavelet"
  return(biwave_object)
}
