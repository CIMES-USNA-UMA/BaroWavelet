#' Continuous Wavelet Transform Transfer Function
#'
#' Computes barorreflex sensitivity as a transfer function suing the Continuous Wavelet Transform
#' @param data A matrix with 3 columns containing time values (first column), RR and SBP values (second
#'             and third column).
#' @param HF Maximum limit of the HF band. Default is 0.4 Hz
#' @param LF Maximum limit of the LF band. Default is 0.15 Hz
#' @param VLF Maximum limit of the VLF band. Default is 0.04 Hz
#' @param chosen.dj Frequency resolution. Default is 1/20
#' @param dt Time resolution (inverse of the sample rate). Default is 0.25
#' @param demean Boolean, should the data be demeaned before analysis? Default is TRUE
#'
#' @return A list with the estimated components of the baroreflex transfer function in the wavelet domain:
#' \item{TransferFun}{The computed baroreflex transfer function}
#' \item{Coherence}{The computed coherence between the two variables}
#' \item{Freqs}{A vector of frequencies for which the transfer function has been computed}
#' \item{Cone}{The computed cone of influence for the Continuous Wavelet Transform}
#' \item{Time}{The original vector of time values}
#' \item{HF}{The chosen maximum limit of the HF band}
#' \item{LF}{The chosen maximum limit of the LF band}
#' \item{VLF}{The chosen maximum limit of the VLF band}
#' \item{type}{A character string specifying which type of transfer function this is}
#'
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @import biwavelet
#' @export
#'
#' @examples
#' # ADD EXAMPLE!
TransferFunCWT <- function(data, HF = 0.4, LF = 0.15, VLF = 0.04,
  chosen.dj = 1/20, dt = 0.25, demean = TRUE){
                  if(demean){
                     for(n in 2:ncol(data)){
                         data[,n] <- data[,n] - mean(data[,n])
                     }
                  }
                  time <- data[,1]
                  WTransform.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                     s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
                  WTransform.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                     s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
                  Smoothed <- SmoothTransforms(WTransform.x, WTransform.y, chosen.dj)
                  sm.WTransform.x <- Smoothed$sm.WTransform.x
                  sm.WTransform.y <- Smoothed$sm.WTransform.y
                  XWTransform <- Smoothed$XWTransform
                  sm.XWTransform <- Smoothed$sm.XWTransform
                  TransferFun <- sm.XWTransform/ sm.WTransform.x
                  Cospectrum <- Re(sm.XWTransform/ sqrt(sm.WTransform.x * sm.WTransform.y))
                  Quadrature <- Im(sm.XWTransform/ sqrt(sm.WTransform.x * sm.WTransform.y))
                  Coherence <- abs(sm.XWTransform)^2 / (sm.WTransform.x * sm.WTransform.y)
                  Phase <- atan2(Quadrature, Cospectrum)
                  return(list(TransferFun = TransferFun, Coherence = Coherence,
                      Freqs = 1/WTransform.x$period, Cone = WTransform.x$coi, Time = data[,1],
                         HF = HF, LF = LF, VLF = VLF, type = "TFun_cwt",
                      Cospectrum = Cospectrum, Quadrature = Quadrature, Phase = Phase,
                      Scales = WTransform.x$scale))
}




#########################################################
# AUXILIARY FUNCTIONS                                   #
#########################################################




SmoothTransforms <- function(x, y, chosen.dj = 1/20){
  N <- nrow(x$wave)
  M <- ncol(x$wave)
  inverse_scales <- matrix(rep(1/t(x$scale), M), ncol = M, 
                           nrow = N)
  XWTransform <- (Conj(x$wave) * y$wave)
  sm.WTransform.x <- biwavelet::smooth.wavelet(inverse_scales * 
                                                 (abs(x$wave)^2), x$dt, chosen.dj, x$scale)
  sm.WTransform.y <- biwavelet::smooth.wavelet(inverse_scales * 
                                                 (abs(y$wave)^2), x$dt, chosen.dj, x$scale)
  sm.XWTransform = biwavelet::smooth.wavelet(inverse_scales * 
                                               XWTransform, x$dt, chosen.dj, x$scale)
  return(list(sm.WTransform.x = sm.WTransform.x, sm.WTransform.y = sm.WTransform.y,
              sm.XWTransform = sm.XWTransform, XWTransform = XWTransform ))
}



GetCWTscales <- function(HF = 0.4, VLF = 0.04, chosen.dj = 1/20){
  max_scale = 1/(VLF - 0.01)
  min_scale = 1/(HF + 0.1)
  limit <- round(max_scale) + 1
  scales = seq(log2(min_scale), log2(limit), by = chosen.dj)
  scales <- scales[scales <= max_scale]
  scales <- 2^scales
  return(scales)
}


GetInverseCWTscales <- function(x, HF = 0.4, VLF = 0.04, chosen.dj = 1/20){
  N <- nrow(x)
  M <- ncol(x)
  scales <- GetCWTscales(HF = HF, VLF = VLF, chosen.dj = chosen.dj)
  inverse_scales <- matrix(rep(1/t(scales), M), ncol = M, 
                           nrow = N)
  return(inverse_scales)
}


