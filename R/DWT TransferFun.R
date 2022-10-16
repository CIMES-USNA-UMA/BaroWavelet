#' Discrete Wavelet Transform Transfer Function
#'
#' Computes barorreflex sensitivity as a transfer function suing the Discrete Wavelet Transform
#' @param data A matrix with 3 columns containing time values (first column), RR and SBP values (second
#'             and third column).
#' @param HF Maximum limit of the HF band. Default is 0.4 Hz
#' @param LF Maximum limit of the LF band. Default is 0.15 Hz
#' @param VLF Maximum limit of the VLF band. Default is 0.04 Hz
#' @param wv A wavelet for the analysis. Default is d4
#' @param hrv Boolean, should the heart rate variability results be returned too? Default is FALSE.
#'
#' @return A list with the estimated components of the baroreflex transfer function in the wavelet domain:
#' \item{HF}{The computed vector of alpha indices for the HF band}
#' \item{LF}{The computed vector of alpha indices for the LF band}
#' \item{VLF}{The chosen maximum limit of the VLF band}
#' \item{Time}{The original vector of time values}
#' \item{type}{A character string specifying which type of transfer function this is}
#' \item{hrv}{A list containing the variability results computed at the HF and LF bands, as well as the
#'            LF/HF ratio. Only returned if the hrv argument is set to TRUE}
#'
#' @details This function works in the following way: first the variables are adapted so that they can be
#' processed by the \href{https://CRAN.R-project.org/package=RHRV}{RHRV} package. Then, it sends the variables to the RHRV functions so that their variability
#' can be calculated using the Pruned-MODWT(PMODWT), an algorithm introduced in package \href{https://CRAN.R-project.org/package=RHRV}{RHRV} package (see the reference
#' section for more details regarding package RHRV and the PMODWT algorithm). The variability results are then used
#' to calculate a time-varying alpha index representative of the baroreflex sensitivity at the HF and LF bands.
#' Also, the function allows to return the computed heart variability results by \href{https://CRAN.R-project.org/package=RHRV}{RHRV} package.
#' @author Alvaro Chao-Ecija
#' @references
#'Leandro Rodriguez-Linares, Xose Vila, Maria Jose Lado, Arturo Mendez, Abraham Otero and Constantino
#'Antonio Garcia (2020). RHRV: Heart Rate Variability Analysis of ECG Data. R package version 4.2.6.
#'https://CRAN.R-project.org/package=RHRV
#'
#'
#'Constantino A. García, Abraham Otero, Xosé Vila, David G. Márquez. A new algorithm for wavelet-based heart rate variability analysis,
#'Biomedical Signal Processing and Control, Volume 8, Issue 6, 2013, Pages 542-550.
#'
#' @import RHRV
#' @export
#'
#' @examples
#' # ADD EXAMPLE!
TransferFunDWT <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
                  wv = "d4", hrv = FALSE){
                  RHRVobjects <- GenerateRHRVObjects(data)
                  RHRVobjects <- SendDataToRHRV(RHRVobjects, HF, LF, VLF)
                  SBP_HF <- RHRVobjects$SBP$FreqAnalysis[[1]]$HF
                  SBP_LF <- RHRVobjects$SBP$FreqAnalysis[[1]]$LF
                  RR_HF <- RHRVobjects$RR$FreqAnalysis[[1]]$HF
                  RR_LF <- RHRVobjects$RR$FreqAnalysis[[1]]$LF
                  output <- list()
                  output$HF <- sqrt(RR_HF / SBP_HF)
                  output$HF[sqrt(SBP_HF) == 0] <- 0
                  output$LF <- sqrt(RR_LF / SBP_LF)
                  output$LF[sqrt(SBP_LF) == 0] <- 0
                  output$Time <- data[,1]
                  output$type <- "TFun_dwt"
                  if(hrv){
                    output$HRV <- list()
                    output$HRV$HF <- RR_HF
                    output$HRV$LF <- RR_LF
                    output$HRV$LFHF <- RR_LF/RR_HF
                  }
                  return(output)
}



#' Generate RHRV objects
#'
#' Organizes the data in an object compatible with package RHRV
#' @param data A matrix with 3 columns containing time values (first column), RR and SBP values (second
#'             and third column).
#'
#' @return A list with two objects, one for each analyzed variable, compatible with the RHRV package functions.
#'
#' @author Alvaro Chao-Ecija
#'
#' @keywords internal
#'
#' @import RHRV
#'
#' @examples
#' # ADD EXAMPLE!
GenerateRHRVObjects <- function(data){
           HRV <- list()
           length(HRV) <- 2
           for(n in 1:2){
               HRV[[n]]$Ext <- "hrv"
               HRV[[n]]$FreqAnalysis <- list()
               HRV[[n]]$Verbose <- FALSE
               HRV[[n]]$Beat <- list()
               HRV[[n]]$Beat$Time <- data[,1]
               HRV[[n]]$Freq_HR <- 1 / abs(diff(data[,1]))[1]
               if(n == 1){
                  HRV[[n]]$HR <- 60000/data[, n + 1]
                  HRV[[n]]$Adapt <- FALSE
               } else {
                  HRV[[n]]$SBP <- data[, n + 1]
                  HRV[[n]]$Adapt <- TRUE
              }
          }
          names(HRV) <- c("RR", "SBP")
          return(HRV)
}

#' Send compatible objects to RHRV
#'
#' Allows compatible objects to be analyzed by the RHRV package.
#' @param RHRVobjects A list of objects generated by \link{WaveletPhysiology}{GenerateRHRVObjects}
#'
#' @return A list with two objects, containing the results from the analyses performed by RHRV.
#'
#' @author Alvaro Chao-Ecija
#'
#' @keywords internal
#'
#' @import RHRV
#'
#' @examples
#' # ADD EXAMPLE!
SendDataToRHRV <- function(RHRVobjects, HF = 0.4, LF = 0.15, VLF = 0.04){
  for(n in 1:2){
    RHRVobjects[[n]] <-
      AdaptSBPtoRHRV(RHRVobjects[[n]])
    RHRVobjects[[n]] <- RHRV::CreateFreqAnalysis(
      RHRVobjects[[n]])
    RHRVobjects[[n]] <- RHRV::CalculatePowerBand(
      RHRVobjects[[n]], ULFmax = 0.01, VLFmin = 0.01,
      VLFmax = VLF, LFmin = VLF,
      LFmax = LF,
      HFmin = LF, HFmax = HF,
      type = "wavelet", wavelet = wv)
  }
  return(RHRVObjects)
}


#' Adapt blood pressure data for RHRV
#'
#' Adapts blood pressure data so that it can be analyzed by the PMODWT.
#' @param RHRVobjects An object to be adapted
#'
#' @return The adapted object
#'
#' @author Alvaro Chao-Ecija
#'
#' @keywords internal
#'
#' @import RHRV
#'
#' @examples
#' # ADD EXAMPLE!
AdaptSBPtoRHRV <- function(RHRVobject){
                  if(RHRVobject$Adapt){
                     RHRVobject$HR <-  60000 / RHRVobject$SBP
                  }
                  return(RHRVobject)
}





