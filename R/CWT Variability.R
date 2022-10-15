VariabilityCWT  <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                                             chosen.dj = 1/20, dt = 0.25, demean = TRUE,
                            smooth = FALSE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  }
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,2]/1000), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  WT.x$power <- s.inv * WT.x$power
  if(smooth){
    WT.x$power <- biwavelet::smooth.wavelet(WT.x$power, 
                                     WT.x$dt, chosen.dj, WT.x$scale) 
  }
  return(list(Variability = WT.x$power.corr, Coherence = NA,
              Freqs = 1/WT.x$period, Cone = WT.x$coi, Time = data[,1],
              HF = HF, LF = LF, VLF = VLF, type = "TFun_cwt", smoothed = smooth))
}

AddCoherenceToVariability <- function(Variability, TFun){
  Variability$Coherence <- TFun$Coherence
  Variability$TransferFun <- Variability$Variability
  return(Variability)
}


PlotCwtVariability <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
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
  plot(fun, plot.sig = FALSE, plot.phase = FALSE, 
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