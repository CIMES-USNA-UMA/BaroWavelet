



#' @export
TestIndRaw <- function(raw, time_flags1, time_flags2, 
                       weight = TRUE, use.thr = TRUE, set.control = 2){
  fun <- list(Time = raw[,"Time"], HR = 60000/raw[,"RR"], SBP = raw[,"SBP"])
  time_flags1 <- time_flags1 * 60
  select_time1 <- fun$Time[(fun$Time >= time_flags1[1]) &
                             (fun$Time <= time_flags1[2])]
  select_time1 <- match(select_time1, fun$Time)
  time_flags2 <- time_flags2 * 60
  select_time2 <- fun$Time[(fun$Time >= time_flags2[1]) &
                             (fun$Time <= time_flags2[2])]
  select_time2 <- match(select_time2, fun$Time)
  if((sum(select_time1 %in% select_time2)) != 0){
    if(((min(select_time1) - min(select_time2)) < 0) & 
       ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[select_time2 < max(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) > 0)){
      select_time2 <- select_time2[select_time2 < min(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[(select_time2 < min(select_time1)) &
                                     (select_time2 > max(select_time1)) ]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      stop("Control interval is inside test interval")
    } else if(((min(select_time1) - min(select_time2)) == 0) & 
              ((max(select_time1) - max(select_time2)) == 0)){
      stop("Intervals are equal")
    }
  }
  deltaHR <- abs(mean(fun$HR[select_time1]) - mean(fun$HR[select_time2]))
  deltaSBP <- abs(mean(fun$SBP[select_time1]) - mean(fun$SBP[select_time2]))
  delta <- c(deltaHR, deltaSBP)
  interval <- c(select_time1, select_time2)
    HR <- fun$HR[interval]
    SBP <- fun$SBP[interval]
  if(set.control == 1){
    Effect <- c(double(NROW(select_time1)), rep(1, NROW(select_time2)))
    Aux <- select_time1
  } else {
    Effect <- c(rep(1, NROW(select_time1)),double(NROW(select_time2)))
    Aux <- select_time2
  }
  HRdata = data.frame(HR = HR, x = fun$Time[interval], Effect = Effect)
  SBPdata = data.frame(SBP = SBP, x = fun$Time[interval], Effect = Effect)
  #HFmodel <- lm(HF ~ 1 + Effect + x + Effect:x, data = HFdata, na.action = na.omit)
  #LFmodel <- lm(LF ~ 1 + Effect + x + Effect:x, data = LFdata, na.action = na.omit)
  HRmodel <- lm(HR ~ 1 + Effect, data = HRdata, na.action = na.omit)
  SBPmodel <- lm(SBP ~ 1 + Effect, data = SBPdata, na.action = na.omit)
  HRtest <- car::linearHypothesis(HRmodel, c("Effect=0"), 
                                  vcov. = sandwich::vcovHAC)
  SBPtest <- car::linearHypothesis(SBPmodel, c("Effect=0"),
                                  vcov. = sandwich::vcovHAC)
  realHR <- lmtest::coeftest(HRmodel, vcov. = sandwich::vcovHAC)
  realSBP <- lmtest::coeftest(SBPmodel, vcov. = sandwich::vcovHAC)
  seHR <- realHR[2,2]
  seSBP <- realSBP[2,2]
  HRp <- HRtest[2,4]
  SBPp <- SBPtest[2,4]
  p.value <- c(HRp, SBPp)
  output <- round(rbind(delta, p.value, c(seHR, seSBP)),4)
  rownames(output)[1] <- "Difference"
  rownames(output)[3] <- "SE"
  colnames(output) <- c("HR", "SBP")
  return(output)
}

#' @export
TestIndHRV <- function(fun, time_flags1, time_flags2, 
                       weight = TRUE, use.thr = TRUE, set.control = 2){
  #if(min(time_flags1) >= max(time_flags2)){
  #  set.control = 1
  #} else if(max(time_flags1) <= min(time_flags2)){
  #  set.control = 2
  #} else {
  #  stop("Overlapping Intervals")
  #}
  time_flags1 <- time_flags1 * 60
  select_time1 <- fun$Time[(fun$Time >= time_flags1[1]) &
                             (fun$Time <= time_flags1[2])]
  select_time1 <- match(select_time1, fun$Time)
  time_flags2 <- time_flags2 * 60
  select_time2 <- fun$Time[(fun$Time >= time_flags2[1]) &
                             (fun$Time <= time_flags2[2])]
  select_time2 <- match(select_time2, fun$Time)
  if((sum(select_time1 %in% select_time2)) != 0){
    if(((min(select_time1) - min(select_time2)) < 0) & 
       ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[select_time2 < max(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) > 0)){
      select_time2 <- select_time2[select_time2 < min(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[(select_time2 < min(select_time1)) &
                                     (select_time2 > max(select_time1)) ]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      stop("Control interval is inside test interval")
    } else if(((min(select_time1) - min(select_time2)) == 0) & 
              ((max(select_time1) - max(select_time2)) == 0)){
      stop("Intervals are equal")
    }
  }
  deltaHF <- abs(mean(fun$HF[select_time1]) - mean(fun$HF[select_time2]))
  deltaLF <- abs(mean(fun$LF[select_time1]) - mean(fun$LF[select_time2]))
  deltaLFHF <- abs(mean(fun$LFHF[select_time1]) - mean(fun$LFHF[select_time2]))
  delta <- c(deltaHF, deltaLF, deltaLFHF)
  interval <- c(select_time1, select_time2)
    HF <- fun$HF[interval]
    LF <- fun$LF[interval]
    LFHF <- fun$LFHF[interval]
    if(set.control == 1){
      Effect <- c(double(NROW(select_time1)), rep(1, NROW(select_time2)))
      Aux <- select_time1
    } else {
      Effect <- c(rep(1, NROW(select_time1)),double(NROW(select_time2)))
      Aux <- select_time2
    }
  HFdata = data.frame(HF = HF, x = fun$Time[interval], Effect = Effect)
  LFdata = data.frame(LF = LF, x = fun$Time[interval], Effect = Effect)
  LFHFdata = data.frame(LFHF = LFHF, x = fun$Time[interval], Effect = Effect)
  #HFmodel <- lm(HF ~ 1 + Effect + x + Effect:x, data = HFdata, na.action = na.omit)
  #LFmodel <- lm(LF ~ 1 + Effect + x + Effect:x, data = LFdata, na.action = na.omit)
  HFmodel <- lm(HF ~ 1 + Effect, data = HFdata, na.action = na.omit)
  LFmodel <- lm(LF ~ 1 + Effect, data = LFdata, na.action = na.omit)
  LFHFmodel <- lm(LFHF ~ 1 + Effect, data = LFHFdata, na.action = na.omit)
  HFtest <- car::linearHypothesis(HFmodel, c("Effect=0"), 
                                  vcov. = sandwich::vcovHAC)
  LFtest <- car::linearHypothesis(LFmodel, c("Effect=0"),
                                  vcov. = sandwich::vcovHAC)
  LFHFtest <- car::linearHypothesis(LFHFmodel, c("Effect=0"),
                                  vcov. = sandwich::vcovHAC)
  realHF <- lmtest::coeftest(HFmodel, vcov. = sandwich::vcovHAC)
  realLF <- lmtest::coeftest(LFmodel, vcov. = sandwich::vcovHAC)
  realLFHF <- lmtest::coeftest(LFHFmodel, vcov. = sandwich::vcovHAC)
  seHF <- realHF[2,2]
  seLF <- realLF[2,2]
  seLFHF <- realLFHF[2,2]
  HFp <- HFtest[2,4]
  LFp <- LFtest[2,4]
  LFHFp <- LFHFtest[2,4]
  p.value <- c(HFp, LFp, LFHFp)
  #output <- round(rbind(abs(BRS1), BRS2, BRS3, p.value),4)
  output <- round(rbind(delta, p.value, c(seHF, seLF, seLFHF)),4)
  #rownames(output)[1:3] <- c("BRS (delta)", "BRS (ratio)", "BRS (%)")
  rownames(output)[1] <- "Difference"
  rownames(output)[3] <- "SE"
  colnames(output) <- c("HF", "LF", "LFHF")
  return(output)
}

SplitByCoherence <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
         weight = FALSE, phase.rest = NULL, use.phase  =FALSE){
  HF <- fun$HF
  LF <- fun$LF
  VLF <- fun$VLF
  fun <- GetBiwaveletObject(fun)
  if(is.null(phase.rest) & use.phase){
    fun$power <- fun$phase
  } 
  if(!is.null(phase.rest)) phase <- fun$phase
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
  if(!is.null(phase.rest)){
    phase.HF <- phase[(freqs <= HF) & (freqs > LF), select_time]
    phase.LF <- phase[(freqs <= LF) & (freqs > VLF), select_time]
  }
  if(use.thr){
    results.HF[which(rsq.HF < thr)] <- NA
    results.LF[which(rsq.LF < thr)] <- NA
  }
  if(!is.null(phase.rest) && phase.rest == "pos"){
    results.HF[which(phase.HF < 0)] <- NA
    results.LF[which(phase.LF < 0)] <- NA
  } else if(!is.null(phase.rest) && phase.rest == "neg"){
    results.HF[which(phase.HF > 0)] <- NA
    results.LF[which(phase.LF > 0)] <- NA
  }
  results.HF <- colMeans(results.HF, na.rm = TRUE)
  results.LF <- colMeans(results.LF, na.rm = TRUE)
  output <- list(HF = results.HF, LF = results.LF)
  return(output)
}

#' @export
TestIndBRS <- function(fun, time_flags1, time_flags2, 
                             weight = TRUE, use.thr = TRUE, set.control = 2){
  Estimate1 <- ExpectedValues(fun, time_flags1, weight = weight)
  Estimate2 <- ExpectedValues(fun, time_flags2, weight = weight)
  BRS1 <- (Estimate2 - Estimate1) #* 100 / Estimate1
  BRS2 <- Estimate2 / Estimate1
  BRS3 <- abs(BRS1) * 100 / Estimate1
  #if(min(time_flags1) >= max(time_flags2)){
  #  set.control = 1
  #} else if(max(time_flags1) <= min(time_flags2)){
  #  set.control = 2
  #} else {
  #  stop("Overlapping Intervals")
  #}
  time_flags1 <- time_flags1 * 60
  select_time1 <- fun$Time[(fun$Time >= time_flags1[1]) &
                             (fun$Time <= time_flags1[2])]
  select_time1 <- match(select_time1, fun$Time)
  time_flags2 <- time_flags2 * 60
  select_time2 <- fun$Time[(fun$Time >= time_flags2[1]) &
                             (fun$Time <= time_flags2[2])]
  select_time2 <- match(select_time2, fun$Time)
  if((sum(select_time1 %in% select_time2)) != 0){
    if(((min(select_time1) - min(select_time2)) < 0) & 
       ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[select_time2 < max(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) > 0)){
      select_time2 <- select_time2[select_time2 < min(select_time1)]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      select_time2 <- select_time2[(select_time2 < min(select_time1)) &
                                     (select_time2 > max(select_time1)) ]
    } else if(((min(select_time1) - min(select_time2)) > 0) & 
              ((max(select_time1) - max(select_time2)) < 0)){
      stop("Control interval is inside test interval")
    } else if(((min(select_time1) - min(select_time2)) == 0) & 
              ((max(select_time1) - max(select_time2)) == 0)){
      stop("Intervals are equal")
    }
  }
  #point <- max(select_time1)
  #if(point == select_time2[1]) select_time2 <- select_time2[-1]
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
  #HFmodel <- lm(HF ~ 1 + Effect + x + Effect:x, data = HFdata, na.action = na.omit)
  #LFmodel <- lm(LF ~ 1 + Effect + x + Effect:x, data = LFdata, na.action = na.omit)
  HFmodel <- lm(HF ~ 1 + Effect, data = HFdata, na.action = na.omit)
  LFmodel <- lm(LF ~ 1 + Effect, data = LFdata, na.action = na.omit)
  HFtest <- car::linearHypothesis(HFmodel, c("Effect=0"), 
                                  vcov. = sandwich::vcovHAC)
  LFtest <- car::linearHypothesis(LFmodel, c("Effect=0"),
                                  vcov. = sandwich::vcovHAC)
  realHF <- lmtest::coeftest(HFmodel, vcov. = sandwich::vcovHAC)
  realLF <- lmtest::coeftest(LFmodel, vcov. = sandwich::vcovHAC)
  seHF <- realHF[2,2]
  seLF <- realLF[2,2]
  HFp <- HFtest[2,4]
  LFp <- LFtest[2,4]
  p.value <- c(HFp, LFp)
  #output <- round(rbind(abs(BRS1), BRS2, BRS3, p.value),4)
  output <- round(rbind(abs(BRS1), p.value, c(seHF, seLF)),4)
  #rownames(output)[1:3] <- c("BRS (delta)", "BRS (ratio)", "BRS (%)")
  rownames(output)[1] <- "BRS (delta)"
  rownames(output)[3] <- "SE"
  return(output)
}