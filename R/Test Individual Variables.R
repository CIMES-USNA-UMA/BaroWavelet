#' Test Individual BRS
#'
#' Significant test to compare two time-domain segments from an individual BRS computation
#' @param fun BRS obtained either by \link[BaroWavelet]{AlphaIndexDWT} or
#' \link[BaroWavelet]{TransferFunCWT}
#' @param time_flags1 A vector containing the minimum and maximum limits for the first time interval, in minutes.
#' @param time_flags2 A vector containing the minimum and maximum limits for the second time interval, in minutes.
#' @param thr Coherence threshold to ensure the reliability of the estimates. Default is 0.5
#' @param use.coherence Boolean, should a coherence threshold be used? Default is TRUE
#' 
#' @return A vector contaning p-values computed for the HF and LF bands
#'
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' AlphaIndex <- AlphaIndexDWT(Data, wv = "d8" error = 0.0005)
#' 
#' TestIndBRS(AlphaIndex, c(0, 1.7), c(35, 38))
#' TestIndBRS(AlphaIndex, c(0, 1.7), c(8, 9.5))
TestIndBRS <- function(fun, time_flags1, time_flags2, thr = 0.5,
                            use.coherence = TRUE){
  if(fun$type == "brs_dwt"){
    pvals <- TestBRSDWT(fun, time_flags1, time_flags2)
  } else if(fun$type == "brs_cwt"){
    pvals <- TestBRSCWT(fun, time_flags1, time_flags2,
                           thr = thr, use.thr = use.coherence)
  }
  return(pvals)
}


TestBRSDWT <- function(fun, time_flags1, time_flags2){
  time_flags1 <- time_flags1 * 60
  time_flags2 <- time_flags2 * 60
  select_time1 <- fun$Time[(fun$Time >= as.integer(time_flags1[1])) &
                              (fun$Time <= as.integer(time_flags1[2]))]
  select_time2 <- fun$Time[(fun$Time >= as.integer(time_flags2[1])) &
                             (fun$Time <= as.integer(time_flags2[2]))]
  select_time1 <- match(select_time1, fun$Time)
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
  HF1 <- fun$HF[select_time1]
  LF1 <- fun$LF[select_time1]
  HF2 <- fun$HF[select_time2]
  LF2 <- fun$LF[select_time2]
  HFtest <- ks.test(HF1, HF2)$p.value
  LFtest <- ks.test(LF1, LF2)$p.value
  output <- c(HFtest, LFtest)
  names(output) <- c("HF", "LF")
  return(output)
}


TestBRSCWT <- function(fun, time_flags1, time_flags2, thr = 0.5, use.thr = TRUE){
  HF <- fun$HF
  LF <- fun$LF
  VLF <- fun$VLF
  time_flags1 <- time_flags1 * 60
  time_flags2 <- time_flags2 * 60
  select_time1 <- fun$Time[(fun$Time >= as.integer(time_flags1[1])) &
                             (fun$Time <= as.integer(time_flags1[2]))]
  select_time2 <- fun$Time[(fun$Time >= as.integer(time_flags2[1])) &
                             (fun$Time <= as.integer(time_flags2[2]))]
  select_time1 <- match(select_time1, fun$Time)
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
  fun <- BaroWavelet:::GetBiwaveletObject(fun)
  freqs <- 1/fun$period
  sel_power <- fun$power
  if(!use.thr) thr <- 0
  results.HF <- fun$power[(freqs <= HF) & (freqs > LF), ]
  results.LF <- fun$power[(freqs <= LF) & (freqs > VLF), ]
  rsq.HF <- fun$rsq[(freqs <= HF) & (freqs > LF), ]
  rsq.LF <- fun$rsq[(freqs <= LF) & (freqs > VLF), ]
  if(use.thr){
    results.HF[which(rsq.HF < thr)] <- NA
    results.LF[which(rsq.LF < thr)] <- NA
  }
  results.HF <- colMeans(results.HF, na.rm = TRUE)
  results.LF <- colMeans(results.LF, na.rm = TRUE)
  HF1 <- results.HF[select_time1]
  LF1 <- results.LF[select_time1]
  HF2 <- results.HF[select_time2]
  LF2 <- results.LF[select_time2]
  HFtest <- ks.test(HF1, HF2)$p.value
  LFtest <- ks.test(LF1, LF2)$p.value
  output <- c(HFtest, LFtest)
  names(output) <- c("HF", "LF")
  return(output)
}


#' Test Individual HR and BP
#'
#' Significant test to compare two time-domain segments in HR and BP recordings
#' @param data A matrix with HR and BP recordings, as well as timestamps
#' @param time_flags1 A vector containing the minimum and maximum limits for the first time interval, in minutes
#' @param time_flags2 A vector containing the minimum and maximum limits for the second time interval, in minutes
#' 
#' @return A vector contaning p-values computed for the HR and the SBP
#'
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' 
#' TestIndBRS(Data, c(0, 1.7), c(35, 38))
#' TestIndBRS(Data, c(0, 1.7), c(8, 9.5))
TestIndHRandBP <- function(data, time_flags1, time_flags2){
  data <- list(Time = data[,"Time"], HR = 60000/data[,"RR"], SBP = data[,"SBP"])
  time_flags1 <- time_flags1 * 60
  select_time1 <- fun$Time[(fun$Time >= as.integer(time_flags1[1])) &
                             (fun$Time <= as.integer(time_flags1[2]))]
  select_time1 <- match(select_time1, fun$Time)
  time_flags2 <- time_flags2 * 60
  select_time2 <- fun$Time[(fun$Time >= as.integer(time_flags2[1])) &
                             (fun$Time <= as.integer(time_flags2[2]))]
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
  HRtest <- ks.test(data$HR[select_time1], data$HR[select_time2])$p.value
  SBPtest <- ks.test(data$SBP[select_time1], data$SBP[select_time2])$p.value
  output <- c(HRtest, SBPtest)
  names(output) <- c("HR", "SBP")
  return(output)
}

#' Test Individual HRV
#'
#' Significant test to compare two time-domain segments from an individual HRV computation
#' @param fun BRS obtainedfrom \link[BaroWavelet]{AlphaIndexDWT}, with HRV data.
#' @param time_flags1 A vector containing the minimum and maximum limits for the first time interval, in minutes.
#' @param time_flags2 A vector containing the minimum and maximum limits for the second time interval, in minutes.
#'
#' @return A vector contaning p-values computed for the HF and LF bands, as well as the LF/HF ratio
#'
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' AlphaIndex <- AlphaIndexDWT(Data, wv = "d8", error = 0.0005, hrv = TRUE)
#' 
#' TestIndHRV(AlphaIndex, c(0, 1.7), c(35, 38))
#' TestIndHRV(AlphaIndex, c(0, 1.7), c(8, 9.5))
TestIndHRV <- function(fun, time_flags1, time_flags2){
  time_flags1 <- time_flags1 * 60
  select_time1 <- fun$Time[(fun$Time >= as.integer(time_flags1[1])) &
                             (fun$Time <= as.integer(time_flags1[2]))]
  select_time1 <- match(select_time1, fun$Time)
  time_flags2 <- time_flags2 * 60
  select_time2 <- fun$Time[(fun$Time >= as.integer(time_flags2[1])) &
                             (fun$Time <= as.integer(time_flags2[2]))]
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
  HFtest <- ks.test(fun$HF[select_time1], fun$HF[select_time2])$p.value
  LFtest <- ks.test(fun$LF[select_time1], fun$LF[select_time2])$p.value
  LFHFtest <- ks.test(fun$LFHF[select_time1], fun$LFHF[select_time2])$p.value
  output <- c(HFtest, LFtest, LFHFtest)
  names(output) <- c("HF", "LF", "LFHF")
  return(output)
}




