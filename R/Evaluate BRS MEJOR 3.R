

EvaluateBRS5 <- function(fun, time_flags1, time_flags2, plot = TRUE, 
   weight = TRUE, use.thr = TRUE, set.control = 1, show.CI = FALSE,
     plotHF = TRUE, plotLF = TRUE, newPlot = TRUE, title){
             Estimate1 <- ExpectedValues(fun, time_flags1, weight = weight)
             Estimate2 <- ExpectedValues(fun, time_flags2, weight = weight)
             BRS1 <- (Estimate2 - Estimate1) #* 100 / Estimate1
             BRS2 <- Estimate2 / Estimate1
             BRS3 <- abs(BRS1) * 100 / Estimate1
             time_flags1 <- time_flags1 * 60
             select_time1 <- fun$Time[(fun$Time >= time_flags1[1]) &
                     (fun$Time <= time_flags1[2])]
             select_time1 <- match(select_time1, fun$Time)
             time_flags2 <- time_flags2 * 60
             select_time2 <- fun$Time[(fun$Time >= time_flags2[1]) &
                     (fun$Time <= time_flags2[2])]
             select_time2 <- match(select_time2, fun$Time)
             point <- max(select_time1)
             if(point == select_time2[1]) select_time2 <- select_time2[-1]
             interval <- c(select_time1, select_time2)
             if(weight){
                     w1 <- 1:NROW(select_time1)
                     w1 <- exp(-((w1 - mean(w1))^2) / (2 * var(w1)))
                     w2 <- 1:NROW(select_time2)
                     w2 <- exp(-((w2 - mean(w2))^2) / (2 * var(w2)))
                     HF <- c(NROW(select_time1) * w1 * fun$HF[select_time1] / 
                        sum(w1), NROW(select_time2) * w2 * fun$HF[select_time2] / 
                           sum(w2))
                     LF <- c(NROW(select_time1) * w1 * fun$LF[select_time1] / 
                        sum(w1), NROW(select_time2) * w2 * fun$LF[select_time2] / 
                           sum(w2))
             } else {
                     HF <- fun$HF[interval]
                     LF <- fun$LF[interval]
             }
             if(set.control == 1){
                Effect <- c(double(NROW(select_time1)), rep(1, NROW(select_time2)))
                Aux <- select_time1
             } else {
                Effect <- c(rep(1, NROW(select_time1)),double(NROW(select_time2)))
                Aux <- select_time2
             }
             HFdata = data.frame(HF = HF, x = fun$Time[interval], Effect = Effect)
             LFdata = data.frame(LF = LF, x = fun$Time[interval], Effect = Effect)
             HFmodel <- lm(HF ~ 1 + Effect, data = HFdata, na.action = na.omit)
             LFmodel <- lm(LF ~ 1 + Effect, data = LFdata, na.action = na.omit)
             HFtest <- car::linearHypothesis(HFmodel, c("Effect=0"), 
               vcov. = sandwich::vcovHAC)
             LFtest <- car::linearHypothesis(LFmodel, c("Effect=0"),
               vcov. = sandwich::vcovHAC)
             realHF <- lmtest::coeftest(HFmodel, vcov. = sandwich::vcovHAC)
             realLF <- lmtest::coeftest(LFmodel, vcov. = sandwich::vcovHAC)
             seHF <- realHF[1,2]
             seLF <- realLF[1,2]
             HFp <- HFtest[2,4]
             LFp <- LFtest[2,4]
             p.value <- c(HFp, LFp)
             if(plot){
                if(newPlot) x11(title = title)
                coef1 <- HFmodel$coefficients[1]
                coef2 <- LFmodel$coefficients[1]
                if(plotHF & plotLF) par(mfrow = c(2,1))
                if(show.CI){
                   ylim1 = c(0, coef1 + 1.96*seHF + 2)
                   ylim2 <- c(0, coef2 + 1.96*seLF + 2)
                } else {
                   ylim1 <- ylim2 <- NULL
                }
                if(plotHF){
                plot(fun$Time, fun$HF, type = "l", ylab = "Baroreflex Activity", 
                     main = "HF", xlab = "Time", ylim = ylim1)               
                lines(fun$Time[select_time1], fun$HF[select_time1], type = "l", ylab = "Baroreflex Activity", 
                     main = "HF", xlab = "Time", col = "blue", ylim = ylim1)
                lines(fun$Time[select_time2], fun$HF[select_time2], type = "l", ylab = "Baroreflex Activity", 
                     main = "HF", xlab = "Time", col = "green", ylim = ylim1)
                lines(fun$Time[interval], predict(HFmodel), col = "red", lwd = 3, ylim = ylim1)
                if(show.CI){
                   lines(fun$Time[interval], rep(coef1 + 1.96*seHF, NROW(interval)),
                      col = "red", ylim = ylim1, lwd = 2)
                   lines(fun$Time[interval], rep(coef1 - 1.96*seHF, NROW(interval)), 
                      col = "red", ylim = ylim1, lwd = 2)
                }
                }
                if(plotLF){
                plot(fun$Time, fun$LF, type = "l", ylab = "Baroreflex Activity", 
                     main = "LF", xlab = "Time", ylim = ylim2)
                lines(fun$Time[select_time1], fun$LF[select_time1], type = "l", ylab = "Baroreflex Activity", 
                     main = "LF", xlab = "Time", col = "blue", ylim = ylim2)
                lines(fun$Time[select_time2], fun$LF[select_time2], type = "l", ylab = "Baroreflex Activity", 
                     main = "LF", xlab = "Time", col = "green", ylim = ylim2)             
                lines(fun$Time[interval], predict(LFmodel), col = "red", lwd = 3, ylim = ylim2)
                if(show.CI){
                   lines(fun$Time[interval], rep(coef2 + 1.96*seLF, NROW(interval)), 
                     col = "red", ylim = ylim2, lwd = 2)
                   lines(fun$Time[interval], rep(coef2 - 1.96*seLF, NROW(interval)), 
                     col = "red", ylim = ylim2, lwd = 2)
                }
                }

             }    
             output <- round(rbind(BRS1, BRS2, BRS3, p.value),4)
             rownames(output)[1:3] <- c("BRS (delta)", "BRS (ratio)", "BRS (%)")
             return(output)
}


ExtractCoherence <- function(fun, time_flags1, time_flags2){
                    Scal1 <- ExtractScalogram(fun, time_flags1)
                    Scal2 <- ExtractScalogram(fun, time_flags2)
                    N1 = ncol(Scal1[[1]])
                    N2 = ncol(Scal2[[1]])
                    if(N1 > N2) for(n in 1:2) Scal1[[n]] <- Scal1[[n]][,(N1-N2+1):N1]
                    if(N2 > N1) for(n in 1:2) Scal2[[n]] <- Scal2[[n]][, 1:N1]
                    coh = list()
                    length(coh) <- 2
                    for(n in 1:2){
                        xx <- Scal1[[n]] * Conj(Scal1[[n]])
                        yy <- Scal2[[n]] * Conj(Scal2[[n]])
                        xy <- Scal1[[n]] * Conj(Scal2[[n]])
                        coh[[n]] = abs(rowMeans(xy))^2 / (rowMeans(xx) * rowMeans(yy))
                    }
                    return(coh)
}


EvaluateBRSinCWT <- function(fun, time_flags1, time_flags2, plot = TRUE, 
   weight = TRUE, use.thr = TRUE){
                  control <- ExtractCWTFreq(fun, time_flags = time_flags1,
                    weight = weight, use.thr = use.thr)
                  test <- ExtractCWTFreq(fun, time_flags = time_flags2,
                    weight = weight, use.thr = use.thr)
                  HF = wilcox.test(control[[1]], test[[1]], paired = TRUE)$p.value
                  LF = wilcox.test(control[[2]], test[[2]], paired = TRUE)$p.value
                  return(round(c(HF, LF),3))
}


ExtractScalogram <- function(fun, time_flags = NULL){
                  HF <- fun$HF
                  LF <- fun$LF
                  VLF <- fun$VLF
                  Time <- fun$Time
                  fun <- GetBiwaveletObject(fun)
                  if(is.null(time_flags)){
                     select_time <- 1:NROW(Time)
                  } else {
                     time_flags <- time_flags * 60
                     select_time <- fun$t[(Time >= time_flags[1]) &
                        (Time <= time_flags[2])]
                     select_time <- match(select_time, Time)
                  }
                  freqs <- 1/fun$period
                  sel_power <- fun$power
                  results.HF <- fun$power[(freqs <= HF) & (freqs > LF), select_time]
                  results.LF <- fun$power[(freqs <= LF) & (freqs > VLF), select_time] 
                  return(list(results.HF, results.LF))
}


ExtractCWTFreq <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
    weight = TRUE){
                  HF <- fun$HF
                  LF <- fun$LF
                  VLF <- fun$VLF
                  fun <- GetBiwaveletObject(fun)
                  if(is.null(time_flags)){
                     select_time <- 1:NROW(fun$t)
                  } else {
                     time_flags <- time_flags * 60
                     select_time <- fun$t[(fun$t >= time_flags[1]) &
                        (fun$t <= time_flags[2])]
                     select_time <- match(select_time, fun$t)
                  }
                  freqs <- 1/fun$period
                  sel_power <- fun$power
                  if(!use.thr) thr <- 0
                  results.HF <- fun$power[(freqs <= HF) & (freqs > LF), select_time]
                  results.LF <- fun$power[(freqs <= LF) & (freqs > VLF), select_time]
                  rsq.HF <- fun$rsq[(freqs <= HF) & (freqs > LF), select_time]
                  rsq.LF <- fun$rsq[(freqs <= LF) & (freqs > VLF), select_time]
                  if(use.thr){
                     results.HF[which(rsq.HF < thr)] <- NA
                     results.LF[which(rsq.LF < thr)] <- NA
                  }        
                  results.HF <- rowMeans(results.HF, na.rm = TRUE)
                  results.LF <- rowMeans(results.LF, na.rm = TRUE)
                  if(weight){
                     tLF <- 1:NROW(results.LF)
                     tHF <- 1:NROW(results.HF)
                     wLF <- exp(-((tLF - mean(tLF))^2) / (2 * var(tLF)))
                     wHF <- exp(-((tHF - mean(tHF))^2) / (2 * var(tHF)))
                  } else {
                     wHF <- rep(1, NROW(results.HF))
                     wLF <- rep(1, NROW(results.LF))
                  }
                  return(list(results.HF*wHF, results.LF*wLF))
}