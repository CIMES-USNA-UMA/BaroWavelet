#' @export
PlotTransferFun <- function(tf, avg = FALSE, time_col = "brown", HFcolor = "yellow",
    LFcolor = "green", time_flags = NULL, thr = 0.5,
       use.thr = TRUE, nfreqs = 7, tem = FALSE, newPlot = TRUE, title, plotHF = TRUE,
    plotLF = TRUE){
        #if(dev.cur() > 1) dev.off()
        if(newPlot) x11(title = paste("Transfer Function from", title))
        if(tf$type == "TFun_dwt"){
           im <- PlotTransferFunDWT(tf, time_flags = time_flags, col = time_col, tem = tem)
        } else if(tf$type == "TFun_cwt"){
           if(avg){
             im <- PlotAvgCwtTransferFun(tf, thr, use.thr, time_flags = time_flags,
              len = nfreqs, HFcolor = HFcolor, LFcolor = LFcolor,
                 Tcolor = time_col, tem = tem)
           } else {
              im <- PlotCwtTransferFun(tf, thr, use.thr, time_flags = time_flags,
                tem = tem)
           }
       }
       if(tem){
          return(im)
       }
}



PlotTransferFunDWT <- function(tf, time_flags = NULL, col = "brown", tem  =FALSE, plotHF = TRUE,
                               plotLF = TRUE, use.xlim = FALSE){
                  time <- tf$Time
                  HF <- tf$HF
                  LF <- tf$LF
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

PlotCwtTransferFun <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
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
                  xlim = NULL
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

PlotAvgCwtTransferFun <- function(fun, thr = 0.5, use.thr = TRUE, scale = 1,
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



AssembleCwtTransferFun <- function(framework, index){
                  Data <- framework$"General Data"
                  TFun <- framework$Analyses[[index]]$BRS$CWT
                  TFun$HF <- Data$HF
                  TFun$LF <- Data$LF
                  TFun$VLF <- Data$VLF
                  TFun$Time <- framework$Analyses[[index]]$Data[,1]
                  return(TFun)
}



