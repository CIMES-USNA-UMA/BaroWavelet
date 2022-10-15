GetCausalCoupling34 <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                        chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5,
                          in_phase = TRUE, out_phase = TRUE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  freqs <- 1/WT.x$period
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  XWT <- (WT.y$wave * Conj(WT.x$wave)) #WT.x is X, WT.y is Y
  # Therefore, argument data is a matrix with three columns that are:
  # [Time, Y, X]. Therefore, if one should study the influence between SBP
  # and RR, then the data should be arranged as [Time, RR, SBP]
  sWT.x <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave)^2), 
                                     WT.x$dt, chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sWT.y <- biwavelet::smooth.wavelet(s.inv * (abs(WT.y$wave)^2), WT.x$dt, 
                                     chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sXWT = biwavelet::smooth.wavelet(s.inv * XWT, WT.x$dt, 
                                   chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence <- abs(sXWT)^2 / (sWT.x * sWT.y)
  Total_GC <- log(1 - Coherence)
  phase <- atan2(Im(sXWT), Re(sXWT))
  phase2 <- phase
  fun11 <- fun12 <- matrix(0, ncol = ncol(phase), nrow = nrow(phase))
  for(n in 1:nrow(phase2)){
    for(m in 1:ncol(phase2)){
      if(((phase[n,m] > 0 & phase[n,m] < pi/2)&in_phase) | ((phase[n,m] > -pi & phase[n,m] < -pi/2)& out_phase)){
        fun11[n,m] <- 1 # X leads, Y lags
      } else if(((phase[n,m] > pi/2 & phase[n,m] < pi)&out_phase) | ((phase[n,m] > -pi/2 & phase[n,m] < 0)& in_phase)){
        fun12[n,m] <- 1 # Y leads, X lags
      }
      if(abs(phase[n,m]) > pi/2){
        phase2[n,m] <- abs(phase2[n,m]) - pi/2
      }
      #if(abs(phase[n,m]) > pi/4){
      #  phase2[n,m] <- abs(phase2[n,m] - pi/2)
      #}
    }
  }
  tol <- tol * pi / 180
  fun2 <- exp((-(phase2 - (pi/4))^2)/(2*tol^2))
  fun2 <- (Coherence >= thr)*1
  sXWT1 = biwavelet::smooth.wavelet(s.inv * XWT * fun11 * fun2, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # X leads, Y lags. X predicts Y
  sXWT2 = biwavelet::smooth.wavelet(s.inv * XWT * fun12 * fun2, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # Y leads, X lags. Y predicts X
  Coherence1 <- abs(sXWT1)^2 / (sWT.x * sWT.y) # X predicts Y
  Coherence2 <- abs(sXWT2)^2 / (sWT.x * sWT.y) # Y predicts X
  Cohesion <- Re(sXWT) / sqrt(sWT.x * sWT.y)
  Cohesion1 <- abs(Re(sXWT1) / sqrt(sWT.x * sWT.y)) # X leads Y (Y causes X)
  Cohesion2 <- abs(Re(sXWT2) / sqrt(sWT.x * sWT.y)) # Y predicts X
  GC_mask1 <- fun11*fun2
  GC_mask2 <- fun12*fun2
  for(n in 1:nrow(Coherence)){
    for(m in 1:ncol(Coherence)){
      if(GC_mask1[n,m] == 0) GC_mask1[n,m] <- 0
      if(GC_mask2[n,m] == 0) GC_mask2[n,m] <- 0
    }
  }
  GC1 <- Total_GC * GC_mask1 # X predicts Y
  GC2 <- Total_GC * GC_mask2 # Y predicts X
  # Remember: X = Output, Y = Input. If you are measuring BRS:
  # X = SBP, Y = RR, and Coherence2 will contain BRS information.
  Objects <- list()
  # The list returned by the funcion will have the following elements: Total Coherence, 
  # Feedback coherence, Feedforward coherence, Total Cohesion, Feedback Cohesion, Feedforward
  # Cohesion
  return(list(Coherence = Coherence, Coherence1 = Coherence1, Coherence2 = Coherence2, Cohesion = Cohesion,
              Cohesion1 = Cohesion1, Cohesion2 = Cohesion2, GC= Total_GC, GC1 = GC1, GC2 = GC2, Freqs = freqs))
}

GetPhaseInducedCoherence <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                                chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5,
                                in_phase = TRUE, out_phase = TRUE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  freqs <- 1/WT.x$period
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  XWT <- (WT.y$wave * Conj(WT.x$wave)) #WT.x is X, WT.y is Y
  # Therefore, argument data is a matrix with three columns that are:
  # [Time, Y, X]. Therefore, if one should study the influence between SBP
  # and RR, then the data should be arranged as [Time, RR, SBP]
  sWT.x <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave)^2), 
                                     WT.x$dt, chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sWT.y <- biwavelet::smooth.wavelet(s.inv * (abs(WT.y$wave)^2), WT.x$dt, 
                                     chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sXWT = biwavelet::smooth.wavelet(s.inv * XWT, WT.x$dt, 
                                   chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence <- abs(sXWT)^2 / (sWT.x * sWT.y)
  Total_GC <- log(1 - Coherence)
  phase <- atan2(Im(sXWT), Re(sXWT))
  phase2 <- phase
  fun11 <- fun12 <- matrix(0, ncol = ncol(phase), nrow = nrow(phase))
  for(n in 1:nrow(phase2)){
    for(m in 1:ncol(phase2)){
      if(((phase[n,m] > 0 & phase[n,m] < pi/2)&in_phase) | ((phase[n,m] > -pi & phase[n,m] < -pi/2)& out_phase)){
        fun11[n,m] <- 1 # X leads, Y lags
      } else if(((phase[n,m] > pi/2 & phase[n,m] < pi)&out_phase) | ((phase[n,m] > -pi/2 & phase[n,m] < 0)& in_phase)){
        fun12[n,m] <- 1 # Y leads, X lags
      }
      if(abs(phase[n,m]) > pi/2){
        phase2[n,m] <- abs(phase2[n,m]) - pi/2
      }
      #if(abs(phase[n,m]) > pi/4){
      #  phase2[n,m] <- abs(phase2[n,m] - pi/2)
      #}
    }
  }
  sXWT1 = biwavelet::smooth.wavelet(s.inv * XWT * fun11, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # X leads, Y lags
  sXWT2 = biwavelet::smooth.wavelet(s.inv * XWT * fun12, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # Y leads, X lags
  Coherence1 <- abs(Im(sXWT1))^2 / (sWT.x * sWT.y) # X leads Y
  Coherence2 <- abs(Im(sXWT2))^2 / (sWT.x * sWT.y) # Y leads X
  Coherence3 <- abs(Re(sXWT))^2 / (sWT.x * sWT.y) # Instantaneous
  # Remember: X = Output, Y = Input. If you are measuring BRS:
  # X = SBP, Y = RR, and Coherence2 will contain BRS information.

  # The list returned by the funcion will have the following elements: Total Coherence, 
  # Feedback coherence, Feedforward coherence, Total Cohesion, Feedback Cohesion, Feedforward
  # Cohesion
  return(list(Coherence = Coherence, Coherence1 = Coherence1, Coherence2 = Coherence2, Coherence3 = Coherence3,
              sXWT = sXWT, sXWT1 = sXWT1,
              sXWT2 = sXWT2, Freqs = freqs))
}



CalculateCoherency <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                                  chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WTC.xy <- biwavelet::wtc(cbind(time, data[,3]), cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01), nrands = 2)
  return(WTC.xy)
}


GetCausalCoupling2 <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                              chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5,
                              in_phase = TRUE, out_phase = TRUE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  freqs <- 1/WT.x$period
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  XWT <- (WT.y$wave * Conj(WT.x$wave)) #WT.x is X, WT.y is Y
  # Therefore, argument data is a matrix with three columns that are:
  # [Time, Y, X]. Therefore, if one should study the influence between SBP
  # and RR, then the data should be arranged as [Time, RR, SBP]
  sWT.x <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave)^2), 
                                     WT.x$dt, chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sWT.y <- biwavelet::smooth.wavelet(s.inv * (abs(WT.y$wave)^2), WT.x$dt, 
                                     chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sXWT = biwavelet::smooth.wavelet(s.inv * XWT, WT.x$dt, 
                                   chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence <- abs(sXWT)^2 / (sWT.x * sWT.y)
  Total_GC <- log(1 - Coherence)
  phase <- atan2(Im(sXWT), Re(sXWT))
  phase2 <- phase
  fun11 <- fun12 <- matrix(0, ncol = ncol(phase), nrow = nrow(phase))
  for(n in 1:nrow(phase2)){
    for(m in 1:ncol(phase2)){
      if(((phase[n,m] > 0 & phase[n,m] < pi/2)&in_phase) | ((phase[n,m] > -pi & phase[n,m] < -pi/2)& out_phase)){
        fun11[n,m] <- 1 # X leads, Y lags
      } else if(((phase[n,m] > pi/2 & phase[n,m] < pi)&out_phase) | ((phase[n,m] > -pi/2 & phase[n,m] < 0)& in_phase)){
        fun12[n,m] <- 1 # Y leads, X lags
      }
      if(abs(phase[n,m]) > pi/2){
        phase2[n,m] <- abs(phase2[n,m]) - pi/2
      }
      #if(abs(phase[n,m]) > pi/4){
      #  phase2[n,m] <- abs(phase2[n,m] - pi/2)
      #}
    }
  }
  tol <- tol * pi / 180
  fun2 <- exp((-(phase2 - (pi/4))^2)/(2*tol^2))
  fun2 <- (Coherence >= thr)*1
  sXWT1 = biwavelet::smooth.wavelet(s.inv * XWT * fun11 * fun2, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # X leads, Y lags. X predicts Y
  sXWT2 = biwavelet::smooth.wavelet(s.inv * XWT * fun12 * fun2, WT.x$dt, 
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  # Y leads, X lags. Y predicts X
  sWT.x2 <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave * fun11*fun2)^2), WT.x$dt, 
                                               chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  coupling <- sWT.x2/sWT.x
  
  return(list(Coupling = coupling, Freqs = freqs))
}
  
GetCausalCoupling4 <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                               chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5,
                               in_phase = TRUE, out_phase = TRUE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                        s0 = 1/HF , max.scale = 1/(VLF ))
  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF ), max.scale = 1/(VLF ))
  freqs <- 1/WT.x$period
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  XWT <- (WT.x$wave * Conj(WT.y$wave)) #WT.x is X, WT.y is Y
  # Therefore, argument data is a matrix with three columns that are:
  # [Time, Y, X]. Therefore, if one should study the influence between SBP
  # and RR, then the data should be arranged as [Time, RR, SBP]
  sWT.x <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave)^2), 
                                     WT.x$dt, chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sWT.y <- biwavelet::smooth.wavelet(s.inv * (abs(WT.y$wave)^2), WT.x$dt, 
                                     chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sXWT = biwavelet::smooth.wavelet(s.inv * XWT, WT.x$dt, 
                                   chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence <- sXWT / sqrt(sWT.x * sWT.y)
  mask = Coherence
  Cospectrum <- Re(Coherence)
  for(n in 1:nrow(Coherence)){
    for(m in 1:ncol(Coherence)){
      if((abs(Coherence[n,m])^2) < thr) mask[n,m] <- NA
    }
  }
  return(list(Cospectrum = Cospectrum, mask = mask, Freqs = freqs))
}

GetCausalCoupling <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                               chosen.dj = 1/20, dt = 0.25, demean = TRUE, tol = 15, thr = 0.5,
                               in_phase = TRUE, out_phase = TRUE){
  if(demean){
    for(n in 2:ncol(data)){
      data[,n] <- data[,n] - mean(data[,n])
    }
  } # Data is a matrix whith data [Time, Output, Input]
  time <- data[,1]
  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                        s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
  freqs <- 1/WT.x$period
  s.inv <- 1/t(WT.x$scale)
  s.inv <- matrix(rep(s.inv, ncol(WT.x$wave)),
                  nrow = nrow(WT.x$wave))
  XWT <- (WT.y$wave * Conj(WT.x$wave)) #WT.x is Input(SBP), WT.y is Output (RR)
  # Thta way we calculate phase(RR) - phase(SBP)
  # Therefore, argument data is a matrix with three columns that are:
  # [Time, Y, X]. Therefore, if one should study the influence between SBP
  # and RR, then the data should be arranged as [Time, RR, SBP]
  sWT.x <- biwavelet::smooth.wavelet(s.inv * (abs(WT.x$wave)^2), 
                                     WT.x$dt, chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sWT.y <- biwavelet::smooth.wavelet(s.inv * (abs(WT.y$wave)^2), WT.x$dt, 
                                     chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  sXWT = biwavelet::smooth.wavelet(s.inv * XWT, WT.x$dt, 
                                   chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence <- abs(sXWT)^2 / (sWT.x * sWT.y)
  phase <- atan2(Im(sXWT), Re(sXWT))
  phase2 <- phase
  fun11 <- fun12 <- matrix(0, ncol = ncol(phase), nrow = nrow(phase))
  for(n in 1:nrow(phase2)){
    for(m in 1:ncol(phase2)){
      if(((phase[n,m] > 0 ))){ # Output leads
        fun11[n,m] <- 1 
      } else if(((phase[n,m] < 0))){  # Input leads
        fun12[n,m] <- 1 
      }
      #if(abs(phase[n,m]) > pi/4){
      #  phase2[n,m] <- abs(phase2[n,m] - pi/2)
      #}
    }
  }
  sXWT1 = biwavelet::smooth.wavelet(s.inv * XWT * fun11, WT.x$dt, # Output leads
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave) 
  sXWT2 = biwavelet::smooth.wavelet(s.inv * XWT * fun12, WT.x$dt, # Input leads
                                    chosen.dj, WT.x$scale) / ncol(WT.x$wave)
  Coherence1 <- abs(Im(sXWT1))^2 / (sWT.x * sWT.y) # Non baroreflex
  Coherence2 <- abs(Im(sXWT2))^2 / (sWT.x * sWT.y) # baroreflex
  Coherence3 <- abs(Re(sXWT2))^2 / (sWT.x * sWT.y) # baroreflex
  mask1 <- mask2 <- matrix(NA, nrow = nrow(Coherence), ncol = ncol(Coherence))
  for(n in 1:nrow(Coherence)){
    for(m in 1:ncol(Coherence)){
      if((Coherence[n,m]) >= thr) mask1[n,m] <- 1
      if((Coherence[n,m]) >= thr) mask2[n,m] <- 1
    }
  }
  # The list returned by the funcion will have the following elements: Total Coherence, 
  # Feedback coherence, Feedforward coherence, Total Cohesion, Feedback Cohesion, Feedforward
  # Cohesion
  return(list(Coherence = Coherence, Coherence1 = Coherence1, Coherence2 = Coherence2, mask1 = mask1, mask2 = mask2,
              Coherence3 = Coherence3,
              Freqs = freqs))
}


DivideByBands <- function(Coupling, HF = 0.4, LF = 0.15, VLF = 0.04, dif = FALSE,
                          coh = TRUE){
  Freqs <- Coupling$Freqs
  Coupling1 <- Coupling[[3]]
  Coupling2 <- Coupling[[2]]
  Coupling3 <- Coupling[[1]]
  HFb <- Freqs[(Freqs <= HF) & (Freqs >= LF)]
  LFb <- Freqs[(Freqs < LF) & (Freqs >= VLF)]
  total <- Freqs[(Freqs <= HF) & (Freqs >= VLF)]
  HF <- match(HFb, Freqs)
  LF <- match(LFb, Freqs)
  total <- match(total, Freqs)
  Coupling1 <- list(Total = colMeans(Coupling1[total,], na.rm = TRUE),
                    HF = colMeans(Coupling1[HF,], na.rm = TRUE), 
                    LF = colMeans(Coupling1[LF,], na.rm = TRUE))
  Coupling2 <- list(Total = colMeans(Coupling2[total,], na.rm = TRUE),
                    HF = colMeans(Coupling2[HF,], na.rm = TRUE), 
                    LF = colMeans(Coupling2[LF,], na.rm = TRUE))
  Coherence <- list(Total = colMeans(Coupling3[total,], na.rm = TRUE),
                    HF = colMeans(Coupling3[HF,], na.rm = TRUE), 
                    LF = colMeans(Coupling3[LF,], na.rm = TRUE))
  Coupling <- list(y.leads.x = Coupling1, x.leads.y = Coupling2, Coherence = Coherence)
  #for(n in 1:3) Coupling[[n]] <- Coupling[[n]] - Coupling2[[n]]*dif
  return(Coupling)
}

