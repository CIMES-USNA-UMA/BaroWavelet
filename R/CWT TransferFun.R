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
                  WT.x <- biwavelet::wt(cbind(time, data[,3]), dj = chosen.dj,
                     s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
                  WT.y <- biwavelet::wt(cbind(time, data[,2]), dj = chosen.dj,
                     s0 = 1/(HF + 0.1), max.scale = 1/(VLF - 0.01))
                  Smoothed <- SmoothTransforms(WT.x, WT.y, chosen.dj)
                  sWT.x <- Smoothed$sWT.x
                  sWT.y <- Smoothed$sWT.y
                  sXWT <- Smoothed$sXWT
                  TransferFun <- sXWT / sWT.x
                  Coherence <- abs(sXWT)^2 / (sWT.x * sWT.y)
                  return(list(TransferFun = TransferFun, Coherence = Coherence,
                      Freqs = 1/WT.x$period, Cone = WT.x$coi, Time = data[,1],
                         HF = HF, LF = LF, VLF = VLF, type = "TFun_cwt"))
}

#########################################################
# AUXILIARY FUNCTIONS                                   #
#########################################################


#' Smooth Continuous Wavelet Transforms results
#'
#' Smooths results obtained for Continuous Wavelet Transform, so as to compute the
#' baroreflex transfer function.
#' @param x The Continuous Wavelet Transform of variable x
#' @param y The Continuous Wavelet Transform of variable y
#' @param chosen.dj Frequency resolution. Default is 1/20

#'
#' @return A list with the smoothed transforms:
#' \item{sWT.x}{The smoothed transform of variable x}
#' \item{sWT.y}{The smoothed transform of variable y}
#' \item{sXWT}{The smoothed cross-wavelet transform between x and y}
#'
#' @details This function is based on the code that function biwavelet uses to calculate the
#' wavelet coherence. For more information, please check the reference section.
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @keywords internal
#'
#'
#' @import biwavelet
#'
#' @examples
#' # ADD EXAMPLE!
SmoothTransforms <- function(x, y, chosen.dj = 1/20){
  inverse_scales <- matrix(rep(1/t(x$scale),ncol(x$wave)), ncol = ncol(x$wave),
                           nrow = nrow(x$wave))
  XWT <- (Conj(x$wave) * y$wave)
  sWT.x <- biwavelet::smooth.wavelet(inverse_scales * (abs(x$wave)^2),
                                     x$dt, chosen.dj, x$scale)
  sWT.y <- biwavelet::smooth.wavelet(inverse_scales * (abs(y$wave)^2), x$dt,
                                     chosen.dj, x$scale)
  sXWT = biwavelet::smooth.wavelet(inverse_scales * XWT, x$dt,
                                   chosen.dj, x$scale)
  return(list(sWT.x = sWT.x, sWT.y = sWT.y, sXWT = sXWT))
}




