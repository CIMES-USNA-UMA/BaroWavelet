#' Calculate individual indices
#'
#' Computes an individual BRS index from a single subject
#' @param fun BRS obtained either by \link[BaroWavelet]{AlphaIndexDWT} or
#' \link[BaroWavelet]{TransferFunCWT}
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL.
#' @param thr Coherence threshold to ensure the reliability of the estimates. Default is 0.5
#' @param use.coherence Boolean, should a coherence threshold be used? Default is TRUE
#' @param method Which method, either the mean or the median, should be used to obtain the indices
#'
#' @return A matrix containing data (either medians and IQRs, or means and standard deviations)
#' identifying indices at the HF and LF bands for the specific subjects.
#'
#' @author Alvaro Chao-Ecija
#'
#' @references
#' Lazaro J, Gil E, Orini M, Laguna P, Bailón R. Baroreflex Sensitivity Measured by Pulse Photoplethysmography. Front Neurosci. 2019 Apr 18;13:339. Available from: https://www.frontiersin.org/articles/10.3389/fnins.2019.00339/full
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' AlphaIndex <- AlphaIndexDWT(Data, wv = "d8", error = 0.0005)
#' 
#' IndividualIndices(AlphaIndex, c(0, 100/60))
IndividualIndices <- function(fun, time_flags = NULL, thr = 0.5,
   use.coherence = TRUE, method = c("median", "mean")){
                  if(fun$type == "brs_dwt"){
                    index <- IndividualIndicesDWT(fun, time_flags, method = method)
                  } else if(fun$type == "brs_cwt"){
                    index <- IndividualIndicesCWT(fun, time_flags = time_flags,
                                           thr = thr, use.thr = use.coherence, method = method)
                  }
                  return(index)
}


IndividualIndicesDWT <- function(fun, time_flags = NULL, method = c("median", "mean")){
  method <- match.arg(method)
  if(method == "mean"){
    method1 <- mean
    method2 <- sd
    output_rnames <- c("mean", "sd")
    
  } else if(method == "median"){
    method1 <- median
    method2 <- quantile
    output_rnames <- c("median", "P25", "P75")
  }
  if(is.null(time_flags)){
    select_time <- 1:NROW(fun$Time)
  } else {
    time_flags <- time_flags * 60
    select_time <- fun$Time[(fun$Time >= time_flags[1]) &
                              (fun$Time <= time_flags[2])]
    select_time <- match(select_time, fun$Time)
  }
  results.HF <- method1(fun$HF[select_time], na.rm = TRUE)
  results.LF <- method1(fun$LF[select_time], na.rm = TRUE)
  results2.HF <- ifelse(method == "mean", method2(fun$HF[select_time], na.rm = TRUE),
                        as.numeric(method2(fun$HF[select_time], na.rm = TRUE, 0.25 )))
  results2.LF <- ifelse(method == "mean", method2(fun$LF[select_time], na.rm = TRUE),
                        as.numeric(method2(fun$LF[select_time], na.rm = TRUE, 0.25 )))
  if(method == "mean"){
     output <- rbind(c(results.HF, results.LF), c(results2.HF, results2.LF))
  } else {
    output <- rbind(c(results.HF, results.LF), c(results2.HF[1], results2.LF[1]), 
                    c(as.numeric(method2(fun$HF[select_time], na.rm = TRUE, 0.75 )), 
                    as.numeric(method2(fun$LF[select_time], na.rm = TRUE, 0.75 ))))
  }
  colnames(output) <- c("HF", "LF")
  rownames(output) <- output_rnames
  return(output)
}


IndividualIndicesCWT <- function(fun, thr = 0.5, use.thr = TRUE, time_flags = NULL,
                               method = c("median", "mean")){
  method <- match.arg(method)
  if(method == "mean"){
    method1 <- mean
    method2 <- sd
    output_rnames <- c("mean", "sd")
    
  } else if(method == "median"){
    method1 <- median
    method2 <- quantile
    output_rnames <- c("median", "P25", "P75")
  }
  HF <- fun$HF
  LF <- fun$LF
  VLF <- fun$VLF
  fun <- BaroWavelet:::GetBiwaveletObject(fun)
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
  results.HF <- method1(results.HF[select_time], na.rm = TRUE)
  results.LF <- method1(results.LF[select_time], na.rm = TRUE)
  results2.HF <- ifelse(method == "mean", method2(results.HF[select_time], na.rm = TRUE),
                        as.numeric(method2(results.HF[select_time], na.rm = TRUE, 0.25 )))
  results2.LF <- ifelse(method == "mean", method2(results.LF[select_time], na.rm = TRUE),
                        as.numeric(method2(results.LF[select_time], na.rm = TRUE, 0.25 )))
  if(method == "mean"){
    output <- rbind(c(results.HF, results.LF), c(results2.HF, results2.LF))
  } else {
    output <- rbind(c(results.HF, results.LF), c(results2.HF[1], results2.LF[1]), 
                    c(as.numeric(method2(results.HF[select_time], na.rm = TRUE, 0.75 )), 
                    as.numeric(method2(results.LF[select_time], na.rm = TRUE, 0.75 ))))
  }
  colnames(output) <- c("HF", "LF")
  rownames(output) <- output_rnames
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
  results.HF <- colMeans(results.HF, na.rm = TRUE)
  results.LF <- colMeans(results.LF, na.rm = TRUE)
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

