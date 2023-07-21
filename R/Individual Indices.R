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
#' @param use.phase Boolean, should indices be computed from CWT phase difference estimates? Default is FALSE
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
IndividualIndices <- function(fun,
                              time_flags = NULL,
                              thr = 0.5,
                              use.coherence = TRUE,
                              method = c("median", "mean"),
                              use.phase = FALSE) {
  if (fun$type == "brs_dwt") {
    index <- IndividualIndicesDWT(fun, time_flags, method = method)
  } else if (fun$type == "brs_cwt") {
    if(!use.phase){
    index <- IndividualIndicesCWT(
      fun,
      time_flags = time_flags,
      thr = thr,
      use.thr = use.coherence,
      method = method
    )
    } else {
      index <- IndividualIndicesPhaseCWT(
        fun,
        time_flags = time_flags,
        thr = thr,
        use.thr = use.coherence,
        method = method
      )
    }
  }
  return(index)
}


IndividualIndicesDWT <-
  function(fun,
           time_flags = NULL,
           method = c("median", "mean")) {
    method <- match.arg(method)
    if (method == "mean") {
      method1 <- mean
      method2 <- sd
      output_rnames <- c("mean", "sd")
      
    } else if (method == "median") {
      method1 <- median
      method2 <- quantile
      output_rnames <- c("median", "P25", "P75")
    }
    if (is.null(time_flags)) {
      select_time <- 1:NROW(fun$Time)
    } else {
      time_flags <- time_flags * 60
      # The following step is done due to floating point precision issues
      limit1 <-
        match(min(abs(time_flags[1] - fun$Time)), abs(time_flags[1] - fun$Time))
      limit2 <-
        match(min(abs(time_flags[2] - fun$Time)), abs(time_flags[2] - fun$Time))
      select_time <- limit1:limit2
    }
    results.HF <- method1(fun$HF[select_time], na.rm = TRUE)
    results.LF <- method1(fun$LF[select_time], na.rm = TRUE)
    results2.HF <-
      ifelse(method == "mean",
             method2(fun$HF[select_time], na.rm = TRUE),
             as.numeric(method2(fun$HF[select_time], na.rm = TRUE, 0.25)))
    results2.LF <-
      ifelse(method == "mean",
             method2(fun$LF[select_time], na.rm = TRUE),
             as.numeric(method2(fun$LF[select_time], na.rm = TRUE, 0.25)))
    if (method == "mean") {
      output <-
        rbind(c(results.HF, results.LF), c(results2.HF, results2.LF))
    } else {
      output <-
        rbind(
          c(results.HF, results.LF),
          c(results2.HF[1], results2.LF[1]),
          c(as.numeric(method2(
            fun$HF[select_time], na.rm = TRUE, 0.75
          )),
          as.numeric(method2(
            fun$LF[select_time], na.rm = TRUE, 0.75
          )))
        )
    }
    colnames(output) <- c("HF", "LF")
    rownames(output) <- output_rnames
    return(output)
  }


IndividualIndicesCWT <-
  function(fun,
           thr = 0.5,
           use.thr = TRUE,
           time_flags = NULL,
           method = c("median", "mean")) {
    method <- match.arg(method)
    if (method == "mean") {
      method1 <- mean
      method2 <- sd
      output_rnames <- c("mean", "sd")
      
    } else if (method == "median") {
      method1 <- median
      method2 <- quantile
      output_rnames <- c("median", "P25", "P75")
    }
    HF <- fun$HF
    LF <- fun$LF
    VLF <- fun$VLF
    fun <- GetBiwaveletObject(fun)
    if (is.null(time_flags)) {
      select_time <- 1:NROW(fun$t)
    } else {
      time_flags <- time_flags * 60
      # The following step is done due to floating point precision issues
      limit1 <-
        match(min(abs(time_flags[1] - fun$t)), abs(time_flags[1] - fun$t))
      limit2 <-
        match(min(abs(time_flags[2] - fun$t)), abs(time_flags[2] - fun$t))
      select_time <- limit1:limit2
    }
    freqs <- 1 / fun$period
    if (!use.thr)
      thr <- 0
    results.HF <- fun$power[(freqs <= HF) & (freqs > LF), ]
    results.LF <- fun$power[(freqs <= LF) & (freqs > VLF), ]
    rsq.HF <- fun$rsq[(freqs <= HF) & (freqs > LF), ]
    rsq.LF <- fun$rsq[(freqs <= LF) & (freqs > VLF), ]
    if (use.thr) {
      results.HF[which(rsq.HF < thr)] <- NA
      results.LF[which(rsq.LF < thr)] <- NA
    }
    results.HF <- colMeans(results.HF, na.rm = TRUE)
    results.LF <- colMeans(results.LF, na.rm = TRUE)
    result.HF <- method1(results.HF[select_time], na.rm = TRUE)
    result.LF <- method1(results.LF[select_time], na.rm = TRUE)
    results2.HF <-
      ifelse(method == "mean",
             method2(results.HF[select_time], na.rm = TRUE),
             as.numeric(method2(results.HF[select_time], na.rm = TRUE, 0.25)))
    results2.LF <-
      ifelse(method == "mean",
             method2(results.LF[select_time], na.rm = TRUE),
             as.numeric(method2(results.LF[select_time], na.rm = TRUE, 0.25)))
    if (method == "mean") {
      output <-
        rbind(c(result.HF, result.LF), c(results2.HF, results2.LF))
    } else {
      output <-
        rbind(
          c(result.HF, result.LF),
          c(results2.HF[1], results2.LF[1]),
          c(as.numeric(
            method2(results.HF[select_time], na.rm = TRUE, 0.75)
          ),
          as.numeric(
            method2(results.LF[select_time], na.rm = TRUE, 0.75)
          ))
        )
    }
    colnames(output) <- c("HF", "LF")
    rownames(output) <- output_rnames
    return(output)
  }



IndividualIndicesPhaseCWT <-
  function(fun,
           thr = 0.5,
           use.thr = TRUE,
           time_flags = NULL,
           method = c("median", "mean")) {
    method <- match.arg(method)
    if (method == "mean") {
      method1 <- mean
      method2 <- sd
      output_rnames <- c("mean", "sd")
      
    } else if (method == "median") {
      method1 <- median
      method2 <- quantile
      output_rnames <- c("median", "P25", "P75")
    }
    HF <- fun$HF
    LF <- fun$LF
    VLF <- fun$VLF
    fun <- GetBiwaveletObject(fun)
    if (is.null(time_flags)) {
      select_time <- 1:NROW(fun$t)
    } else {
      time_flags <- time_flags * 60
      # The following step is done due to floating point precision issues
      limit1 <-
        match(min(abs(time_flags[1] - fun$t)), abs(time_flags[1] - fun$t))
      limit2 <-
        match(min(abs(time_flags[2] - fun$t)), abs(time_flags[2] - fun$t))
      select_time <- limit1:limit2
    }
    freqs <- 1 / fun$period
    if (!use.thr)
      thr <- 0
    results.HF <- fun$phase[(freqs <= HF) & (freqs > LF), ]
    results.LF <- fun$phase[(freqs <= LF) & (freqs > VLF), ]
    rsq.HF <- fun$rsq[(freqs <= HF) & (freqs > LF), ]
    rsq.LF <- fun$rsq[(freqs <= LF) & (freqs > VLF), ]
    if (use.thr) {
      results.HF[which(rsq.HF < thr)] <- NA
      results.LF[which(rsq.LF < thr)] <- NA
    }
    results.HF <- colMeans(results.HF, na.rm = TRUE)
    results.LF <- colMeans(results.LF, na.rm = TRUE)
    result.HF <- method1(results.HF[select_time], na.rm = TRUE)
    result.LF <- method1(results.LF[select_time], na.rm = TRUE)
    results2.HF <-
      ifelse(method == "mean",
             method2(results.HF[select_time], na.rm = TRUE),
             as.numeric(method2(results.HF[select_time], na.rm = TRUE, 0.25)))
    results2.LF <-
      ifelse(method == "mean",
             method2(results.LF[select_time], na.rm = TRUE),
             as.numeric(method2(results.LF[select_time], na.rm = TRUE, 0.25)))
    if (method == "mean") {
      output <-
        rbind(c(result.HF, result.LF), c(results2.HF, results2.LF))
    } else {
      output <-
        rbind(
          c(result.HF, result.LF),
          c(results2.HF[1], results2.LF[1]),
          c(as.numeric(
            method2(results.HF[select_time], na.rm = TRUE, 0.75)
          ),
          as.numeric(
            method2(results.LF[select_time], na.rm = TRUE, 0.75)
          ))
        )
    }
    colnames(output) <- c("HF", "LF")
    rownames(output) <- output_rnames
    return(output)
  }

#' Calculate time domain values
#'
#' Computes HR and BP time domain indices
#' @param data matrix with 3 columns containing time values (first column), RR and SBP values (second
#'             and third column), or a data frame containing these variables.
#' @param time_flags A vector containing the minimum and maximum limits of a time interval, in minutes.
#'                   Default is NULL.
#' @param method Which method, either the mean or the median, should be used to obtain the indices
#'
#' @return A list with HR and BP levels
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#'
#' TimeDomainValues(Data, c(0, 100/60), "mean")
TimeDomainValues <-
  function(data,
           time_flags = NULL,
           method = c("median", "mean")) {
    if (!is.data.frame(data))
      data <- as.data.frame(data)
    method <- match.arg(method)
    if (is.null(time_flags)) {
      select_time <- 1:NROW(data$Time)
    } else {
      time_flags <- time_flags * 60
      # The following step is done due to floating point precision issues
      limit1 <-
        match(min(abs(time_flags[1] - data$Time)), abs(time_flags[1] - data$Time))
      limit2 <-
        match(min(abs(time_flags[2] - data$Time)), abs(time_flags[2] - data$Time))
      select_time <- limit1:limit2
    }
    if (method == "mean") {
      HR <- SBP <- double(2)
      names(HR) <- names(SBP) <- c("Mean", "SD")
      HR[1] <- mean(60000 / data$RR[select_time])
      HR[2] <- sd(60000 / data$RR[select_time])
      SBP[1] <- mean(data$SBP[select_time])
      SBP[2] <- sd(data$SBP[select_time])
    } else {
      HR <- SBP <- double(3)
      names(HR) <-
        names(SBP) <- c("Median", "P25", "P75")
      HR[1] <- median(60000 / data$RR[select_time])
      HR[2] <-
        as.numeric(quantile(60000 / data$RR[select_time], 0.25))
      HR[3] <-
        as.numeric(quantile(60000 / data$RR[select_time], 0.75))
      SBP[1] <- median(data$SBP[select_time])
      SBP[2] <-
        as.numeric(quantile(data$SBP[select_time], 0.25))
      SBP[3] <-
        as.numeric(quantile(data$SBP[select_time], 0.75))
    }
    return(list(HR = HR, SBP = SBP))
  }

