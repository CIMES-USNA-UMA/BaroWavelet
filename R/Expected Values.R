#' Calculate expected values
#'
#' Computes the expected values of a transfer function
#' @param fun A transfer function obtained either by \link[WaveletPhysiology]{TransferFunDWT} or
#' \link[WaveletPhysiology]{TransferFunCWT}
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL.
#' @param weight Boolean, should results be weighted by a Gaussian function? Default is FALSE.
#' @param thr Coherence threshold to ensure the reliability of the estimates. Default is 0.5
#' @param use.coherence Boolean, should a coherence threshold be used? Default is TRUE
#'
#' @return A vector containing the computed expected values at the HF and LF bands.
#'
#' @author Alvaro Chao-Ecija
#'
#' @references
#'Violeta I. McLoone, John Ringwood. A system identification approach to baroreflex sensitivity
#'estimation. ISSC 2012; 2012.p. 1-6.
#'
#' @export
#'
#' @examples
#' # ADD EXAMPLE!
ExpectedValues <- function(fun, time_flags = NULL, weight = TRUE, thr = 0.5,
   use.coherence = TRUE){
                  if(fun$type == "TFun_dwt"){
                    eVals <- ExpectedValsDWT(fun, time_flags, weight)
                  } else if(fun$type == "TFun_cwt"){
                    eVals <- ExpectedValsCWT(fun, time_flags = time_flags,
                     weight = weight, thr = thr, use.thr = use.coherence)
                  }
                  return(eVals)
}

ExpectedValsDWT <- function(fun, time_flags = NULL, weight = TRUE){
                  if(is.null(time_flags)){
                     select_time <- 1:NROW(fun$Time)
                  } else {
                     time_flags <- time_flags * 60
                     select_time <- fun$Time[(fun$Time >= time_flags[1]) &
                        (fun$Time <= time_flags[2])]
                     select_time <- match(select_time, fun$Time)
                  }
                  if(weight){
                     w <- 1:NROW(select_time)
                     w <- exp(-((w - mean(w))^2) / (2 * var(w)))
                  } else {
                     w <- rep(1,NROW(select_time))
                  }
                  results.HF <- weighted.mean(fun$HF[select_time], w, na.rm = TRUE)
                  results.LF <- weighted.mean(fun$LF[select_time], w, na.rm = TRUE)
                  output <- c(results.HF, results.LF)
                  names(output) <- c("HF", "LF")
                  return(output)
}



ExpectedValsCWT <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
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
                  results.HF <- weighted.mean(results.HF, wHF, na.rm = TRUE)
                  results.LF <- weighted.mean(results.LF, wLF, na.rm = TRUE)
                  output <- c(results.HF, results.LF)
                  names(output) <- c("HF", "LF")
                  return(output)
}

ExpectedPhaseCWT <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
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
  sel_phase <- fun$phase
  if(!use.thr) thr <- 0
  results.HF <- fun$phase[(freqs <= HF) & (freqs > LF), select_time]
  results.LF <- fun$phase[(freqs <= LF) & (freqs > VLF), select_time]
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
  results.HF <- weighted.mean(results.HF, wHF, na.rm = TRUE)
  results.LF <- weighted.mean(results.LF, wLF, na.rm = TRUE)
  output <- c(results.HF, results.LF)
  names(output) <- c("HF", "LF")
  return(output)
}


TimeDomainValues <- function(data, time_flags = NULL){
                  if(is.null(time_flags)){
                     select_time <- 1:NROW(data[,1])
                  } else {
                     time_flags <- time_flags * 60
                     select_time <- data[,1][(data[,1] >= time_flags[1]) &
                        (data[,1] <= time_flags[2])]
                     select_time <- match(select_time, data[,1])
                  }
                  HR <- SBP <- double(2)
                  names(HR) <- names(SBP) <- c("Mean", "SD")
                  HR[1] <- mean(60000/data[,2][select_time])
                  HR[2] <- sd(60000/data[,2][select_time])
                  SBP[1] <- mean(data[,3][select_time])
                  SBP[2] <- sd(data[,3][select_time])
                  return(list(HR = HR, SBP = SBP))
}

