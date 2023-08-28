

#' Plot individual HRV
#'
#' Plots HRV data computed by function \link[RHRV]{CalculatePowerBand} from package \href{https://CRAN.R-project.org/package=RHRV}{RHRV}
#' @param hrv HRV obtained either by RHRV function \link[RHRV]{CalculatePowerBand} or from \link[BaroWavelet]{AlphaIndexDWT}. If
#'            the former method is used, it must be a list with variables named HF, LF and LFHF. If the latter method is
#'            used, argument hrv must be set to TRUE.
#' @param time Time values for the HRV data
#' @param col Color used to highlight a specific time interval. Default is brown
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL.
#' @param use.xlim Boolean. Should the argument time_flags be used as limits of the x-axis? Default is FALSE.
#' @param tem Boolean, creates a temporal file for the plot. Default is FALSE
#' @param newPlot Boolean, generates a new plot without overwriting a previous plot. Default is TRUE
#' @param title Part of the title for the plot. Used to indicate the data source on the title. Default is "data"
#' @param plotHF Boolean, plot results form the HF band. Default is TRUE
#' @param plotLF Boolean, plot results from the LF band. Default is TRUE
#' @param ratio Boolean. Should the LF/HF ratio be plotted? Default is FALSE. Arguments plotHF and plotLF must also be
#'              set to TRUE.
#' @param ylim Maximum y axis limit. Default is NULL
#' @param use.ggplot Boolean, use methods from \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2} package to plot the results. Default is TRUE
#' @param fHF Maximum limit of the HF band, shown at the plot title. Default is 0.4 Hz
#' @param fLF Maximum limit of the LF band, shown at the plot title. Default is 0.15 Hz
#' @param fVLF Maximum limit of the VLF band, shown at the plot title. Default is 0.04 Hz
#' @param size.axis Percentage of scaling of axis values. Default is 100
#' @param size.labels Percentage of scaling of axis labels. Default is 100
#' @param size.title Percentage of scaling of plot titles. Default is 100
#' @param d.labels Percentage of scaling of distance between axis and labels. Default is 100
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
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' AlphaIndex <- AlphaIndexDWT(Data, wv = "d8", error = 0.0005, hrv = TRUE)
#'
#'
#' PlotHRV(AlphaIndex$HRV, AlphaIndex$Time, plotHF = FALSE, plotLF = TRUE,
#' newPlot = FALSE)
#'

PlotHRV <-
  function(hrv,
           time,
           time_flags = NULL,
           col = "brown",
           use.xlim = FALSE,
           tem  = FALSE,
           newPlot = TRUE,
           title = "data",
           plotHF = TRUE,
           plotLF = TRUE,
           ratio = TRUE,
           ylim = NULL,
           use.ggplot = TRUE,
           fHF = 0.4,
           fLF = 0.15,
           fVLF = 0.04,
           size.axis = 100,
           size.labels = 100,
           size.title = 100,
           d.labels = 100) {
    if (newPlot & !tem) {
      dev.new(title = paste("Heart Rate Variability from", title))
    } else if (!newPlot & !tem & !use.ggplot) {
      if(!is.null(dev.list())) dev.off()
    }
    HF <- hrv$HF
    LF <- hrv$LF
    LFHF <- hrv$LFHF
    size.axis <- size.axis / 100
    size.labels <- size.labels / 100
    size.title <- size.title / 100
    d.labels <- d.labels * 3 / 100
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
      plottingData <-
        data.frame(
          time = time,
          LF = LF,
          HF = HF,
          LFHF = LFHF
        )
      if (plotHF & !plotLF) {
        im <-
          ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = time, y = HF)) +
          ggplot2::geom_line()
        if (use.xlim)
          im <-
            im + ggplot2::xlim(time_flags[1] * 60, time_flags[2] * 60)
        return(im)
      } else if (plotLF & !plotHF) {
        im <-
          ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = time, y = LF)) +
          ggplot2::geom_line()
        if (use.xlim)
          im <-
            im + ggplot2::xlim(time_flags[1] * 60, time_flags[2] * 60)
        return(im)
      } else if (plotLF & plotHF & ratio) {
        im <-
          ggplot2::ggplot(data = plottingData,
                          mapping = ggplot2::aes(x = time, y = LFHF)) +
          ggplot2::geom_line()
        if (use.xlim)
          im <-
            im + ggplot2::xlim(time_flags[1] * 60, time_flags[2] * 60)
        return(im)
      } else if (plotLF & plotHF) {
        stop("Not yet implemented")
      }
      
      if (tem) {
        dev.off()
        return(im)
      }
    } else {
      if (plotHF & plotLF & !ratio)
        par(mfrow = c(1, 2), mgp = c(d.labels, 1, 0))
      if (((plotHF & !plotLF) |  (!plotHF & plotLF)) && !ratio)
        par(mgp = c(d.labels, 1, 0))
      if (plotHF & plotLF & ratio) {
        par(mgp = c(d.labels, 1, 0))
        plot(
          time,
          LFHF,
          type = "l",
          xlab = "time (s)",
          ylab = "LF/HF ratio",
          main = "HRV: LF/HF ratio",
          ylim = if (!is.null(ylim))
            c(0, ylim),
          cex.axis = size.axis,
          cex.lab = size.labels,
          cex.main = size.title
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
            band <- LFHF
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
      } else {
        if (plotHF) {
          plot(
            time,
            HF,
            type = "l",
            xlab = "time (s)",
            ylab = expression(HRV ~ ms ^ 2),
            main = paste("HF band (", fHF, " - ", fLF, " Hz)", sep = ""),
            ylim = if (!is.null(ylim))
              c(0, ylim),
            cex.axis = size.axis,
            cex.lab = size.labels,
            cex.main = size.title
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
            xlab = "time (s)",
            ylab = expression(HRV ~ ms ^ 2),
            main = paste("LF band (", fLF, " - ", fVLF, " Hz)", sep = ""),
            ylim = if (!is.null(ylim))
              c(0, ylim),
            cex.axis = size.axis,
            cex.lab = size.labels,
            cex.main = size.title
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
    
  }
