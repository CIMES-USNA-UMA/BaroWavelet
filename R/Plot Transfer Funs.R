#' Plot individual BRS
#'
#' Plots an already computed BRS from a specific subject
#' @param fun BRS obtained either by \link[BaroWavelet]{AlphaIndexDWT} or
#' \link[BaroWavelet]{TransferFunCWT}
#' @param avg Boolean. Plot the scale-averaged CWT BRS? Default is FALSE
#' @param time_col Color used to highlight a specific time interval. Default is brown
#' @param HFcolor Color to be used to highlight the HF band
#' @param LFcolor Color to be used to highlight the LF band
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL.
#' @param thr Coherence threshold to be used for the plot. Default is 0.5
#' @param use.thr Boolean, should a coherence threshold be used in the analyses? Default is TRUE.
#' @param nfreqs Number of frequencies to be used. Default is 7
#' @param tem Boolean, creates a temporal file for the plot. Default is FALSE
#' @param newPlot Boolean, generates a new plot without overwriting a previous plot. Default is TRUE
#' @param title Part of the title for the plot. Used to indicate the data source on the title. Default is "data"
#' @param plotHF Boolean, plot results form the HF band. Default is TRUE
#' @param plotLF Boolean, plot results from the LF band. Default is TRUE
#'
#' @return None
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' Study <- BuildStructure()
#' Study <- AddAnalysis(Study, name = "Simulation")
#' Study <- AddDataToAnalysis(Study, 1, Data$RR, Data$SBP, Data$Time)
#' Study <- AnalyzeBRS(Study, 1)
#' Study <- GetAvgCwtData(Study, 1)
#' 
#' 
#' PlotAnalyzedBRS(Study, 1, "dwt")
#' 
#' PlotAnalyzedBRS(Study, 1, "cwt")
PlotBRS <- function(fun, avg = FALSE, time_col = "brown", HFcolor = "yellow",
    LFcolor = "green", time_flags = NULL, thr = 0.5,
       use.thr = TRUE, nfreqs = 7, tem = FALSE, newPlot = TRUE, title = "data", plotHF = TRUE,
    plotLF = TRUE){
        #if(dev.cur() > 1) dev.off()
        if(newPlot) x11(title = paste("Transfer Function from", title))
        if(fun$type == "brs_dwt"){
           im <- PlotDwtBRS(fun, time_flags = time_flags, col = time_col, tem = tem)
        } else if(fun$type == "brs_cwt"){
           if(avg){
             im <- PlotAvgCwtBRS(fun, thr, use.thr, time_flags = time_flags,
              len = nfreqs, HFcolor = HFcolor, LFcolor = LFcolor,
                 Tcolor = time_col, tem = tem)
           } else {
              im <- PlotCwtBRS(fun, thr, use.thr, time_flags = time_flags,
                tem = tem)
           }
       }
       if(tem){
          return(im)
       }
}



PlotDwtBRS <- function(fun, time_flags = NULL, col = "brown", tem  =FALSE, plotHF = TRUE,
                               plotLF = TRUE, use.xlim = FALSE){
                  time <- fun$Time
                  HF <- fun$HF
                  LF <- fun$LF
                  if(tem){
                     im <- tempfile(fileext = ".png")
                     png(filename = im, width = 6, height = 6, units = "in", res = 400)
                  }
                  plottingData <- data.frame(Time = time, LF = LF, HF = HF)
                  if(plotHF & plotLF) par(mfrow = c(2,1))
                  if(plotHF){
                  im <- ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = HF)) +
                      ggplot2::geom_line()
                  if(use.xlim) im <- im + ggplot2::xlim(time_flags[1]*60, time_flags[2]*60)
                  return(im)
                  }else{
                  im <- ggplot2::ggplot(data = plottingData, mapping = ggplot2::aes(x = Time, y = LF)) +
                      ggplot2::geom_line()
                  if(use.xlim) im <- im + ggplot2::xlim(time_flags[1]*60, time_flags[2]*60)
                  return(im)
                  }
                  if(tem){
                     dev.off()
                     return(im)
                  }

}

PlotCwtBRS <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
   tem  = FALSE, Max = NULL, use.xlim = FALSE, show.coi = TRUE){
                  HF <- fun$HF
                  LF <- fun$LF
                  VLF <- fun$VLF
                  if(use.thr){
                     mask <- plot.contour <- TRUE
                  } else {
                     mask <- plot.contour <- FALSE
                     thr <- 0
                  }
                  if(is.null(Max)) Max <- max(abs(fun$TransferFun))
                  fun <- GetBiwaveletObject(fun, thr = thr)
                  if(mask){
                     fun$wave[which(fun$rsq < thr)] <- -0.1
                  }
                  if(plot.contour){
                     col <- c("white", "#00007F", "blue", "#007FFF",
                        "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
                  } else {
                     col <- NULL
                  }
                  if(tem){
                     im <- tempfile(fileext = ".png")
                     png(filename = im, width = 1500, height = 400)
                  }
                  xlim <- NULL
                  if(use.xlim) xlim = time_flags * 60
                  par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
                  biwavelet::plot.biwavelet(fun, plot.sig = FALSE, plot.phase = TRUE,
                     type = "wavelet", plot.cb = TRUE, arrow.cutoff = thr, plot.coi = show.coi,
                        main = "Transfer Function by CWT (ms/mmHg)",
                          ylim = c(1/HF, 1/VLF), xlim = xlim, fill.cols = col, zlim = c(0,Max))
                  if(plot.contour){
                    contour(fun$t, log2(fun$period), t(fun$rsq), level = thr,
                       col = "black", lwd = 5, add = TRUE, drawlabels = FALSE)
                  }
                  if(!is.null(time_flags)){
                     time_flags = time_flags * 60
                     abline(v = time_flags[1], lwd = 3, col = "black")
                     abline(v = time_flags[2], lwd = 3, col = "black")
                  }
                  abline(h = log2(1/LF), lwd = 3, col = "brown")
                  if(tem){
                     dev.off()
                     return(im)
                  }
}

PlotAvgCwtBRS <- function(fun, thr = 0.5, use.thr = TRUE, scale = 1,
    time_flags = NULL, len = 7, HFcolor = "yellow", LFcolor = "green",
      Tcolor = "brown", tem  =FALSE){
                  HF <- fun$HF
                  LF <- fun$LF
                  VLF <- fun$VLF
                  fun <- GetBiwaveletObject(fun)
                  freqs <- 1/fun$period
                  sel_power <- fun$power
                  if(use.thr){
                     fun$power[which(fun$rsq < thr)] <- NA
                  }
                  min_freqs <- min(freqs)
                  max_freqs <- max(freqs)
                  if(max_freqs > HF) max_freqs <- HF
                  if(min_freqs < VLF) min_freqs <- VLF
                  if(tem){
                     im <- tempfile(fileext = ".png")
                     png(filename = im, width = 6, height = 6, units = "in", res = 400)
                  }
                  par(mfrow = c(3,1))
                  if(!is.null(time_flags)){
                     time_flags <- time_flags * 60
                     select_t <- fun$t[(fun$t >= time_flags[1]) & (fun$t <= time_flags[2])]
                     freq_results <- rowMeans(fun$power[,match(select_t, fun$t)], na.rm = TRUE)
                  } else {
                     freq_results <- rowMeans(fun$power, na.rm = TRUE)
                  }
                  time_results.HF <- colMeans(fun$power[(freqs <= HF) & (freqs > LF),],
                     na.rm = TRUE)
                  time_results.LF <- colMeans(fun$power[(freqs <= LF) & (freqs > VLF),],
                     na.rm = TRUE)
                  freq_results[is.na(freq_results)] <- 0
                  time_results.HF[is.na(time_results.HF)] <- 0
                  time_results.LF[is.na(time_results.LF)] <- 0
                  plot(log2(1/fun$period), freq_results,
                     "l", xlab = "Frequency", ylab = "BRS (ms/mmHg)",
                        main = "Frequency Domain Transfer Function", ylim =
                          c(0, max(freq_results)), xaxt = "n")
                  polygon(x = c(log2(freqs[(freqs <= HF) & (freqs >= VLF)]),
                      rev( log2(freqs[(freqs <= HF) & (freqs >= VLF)]))),
                    y = c(freq_results[(freqs <= HF) & (freqs >= VLF)],
                      double(NROW(freqs[(freqs <= HF) & (freqs >= VLF)]))),
                      col = HFcolor)
                  polygon(x = c(log2(freqs[(freqs <= LF) & (freqs > VLF)]),
                      rev( log2(freqs[(freqs <= LF) & (freqs > VLF)]))),
                    y = c(freq_results[(freqs <= LF) & (freqs > VLF)],
                      double(NROW(freqs[(freqs <= LF) & (freqs > VLF)]))),
                      col = LFcolor)
                  axis(1, at = log2(1/fun$period[seq(1,NROW(fun$period), len = len)]),
                     labels = round(1/fun$period[seq(1,NROW(fun$period), len = len)],3))
                  plot(fun$t, time_results.HF, "l", xlab = "Time", ylab = "BRS (ms/mmHg)",
                        main = "Time Domain Transfer Function (HF band)")
                  if(!is.null(time_flags)){
                     polygon(x = c(select_t, rev(select_t)),
                       y = c(time_results.HF[match(select_t, fun$t)],
                      double(NROW(select_t))),
                      col = Tcolor)
                  }
                  plot(fun$t, time_results.LF, "l", xlab = "Time", ylab = "BRS (ms/mmHg)",
                        main = "Time Domain Transfer Function (LF band)")
                  if(!is.null(time_flags)){
                     polygon(x = c(select_t, rev(select_t)),
                       y = c(time_results.LF[match(select_t, fun$t)],
                      double(NROW(select_t))),
                      col = Tcolor)
                  }
                  if(tem){
                     dev.off()
                     return(im)
                  }

}

GetBiwaveletObject <- function(data, use.thr = TRUE, thr = 0.5){
                  biwave_object <- list()
                  biwave_object$type <- "wtc"
                  biwave_object$t <- data$Time
                  biwave_object$period <- 1/data$Freqs
                  biwave_object$power <- abs(data$TransferFun)
                  biwave_object$wave <- abs(data$TransferFun)
                  #biwave_object$phase <- atan2(Im(data$TransferFun),
                  #   Re(data$TransferFun))
                  biwave_object$phase <- data$Phase
                  biwave_object$rsq <- data$Coherence
                  biwave_object$coi <- data$Cone
                  if(use.thr){
                     biwave_object$signif <- 5 * (biwave_object$rsq >= thr)
                  }
                  class(biwave_object) <- "biwavelet"
                  return(biwave_object)
}



AssembleCwtBRS <- function(framework, index){
                  Data <- framework$"General Data"
                  funun <- framework$Analyses[[index]]$BRS$CWT
                  funun$HF <- Data$HF
                  funun$LF <- Data$LF
                  funun$VLF <- Data$VLF
                  funun$Time <- framework$Analyses[[index]]$Data[,1]
                  return(funun)
}



