#' Continuous Wavelet Transform Transfer Function
#'
#' Computes barorreflex sensitivity as a transfer function using the Continuous Wavelet Transform
#' @param data A matrix with 3 columns containing time values (first column), RR and SBP values (second
#'             and third column).
#' @param HF Maximum limit of the HF band. Default is 0.4 Hz
#' @param LF Maximum limit of the LF band. Default is 0.15 Hz
#' @param VLF Maximum limit of the VLF band. Default is 0.04 Hz
#' @param chosen.dj Frequency resolution to be passed to the \code{link[biwavelet]{wt}} function. Default is 1/20
#' @param demean Boolean, should the data be demeaned before analysis? Default is TRUE
#' @param smooth Boolean, smooth the transforms so that the BRS correlates with the coherence. Default is TRUE
#' @param alpha Boolean, compute an alpha index representative of the BRS instead of the transfer function. Default is FALSE
#' @return A list with the estimated components of the baroreflex transfer function in the wavelet domain:
#' \item{TransferFun}{The computed baroreflex transfer function}
#' \item{Coherence}{The computed coherence between the two variables}
#' \item{Freqs}{A vector of frequencies for which the transfer function has been computed}
#' \item{Cone}{The computed cone of influence for the Continuous Wavelet Transform}
#' \item{Time}{The original vector of time values}
#' \item{HF}{The chosen maximum limit of the HF band}
#' \item{LF}{The chosen maximum limit of the LF band}
#' \item{VLF}{The chosen maximum limit of the VLF band}
#' \item{type}{A character string specifying which type of BRS has been computed}
#'
#' @details This function makes use of the Continuous Wavelet Transform methods
#' provided by package \href{https://CRAN.R-project.org/package=biwavelet}{biwavelet} to
#' compute the baroreflex sensitivity. It employs the functions \code{link[biwavelet]{wt}}
#' and \code{link[biwavelet]{smooth.wavelet}}. This last function is used in a smoothing routine
#' based on the one that biwavelet function \code{link[biwavelet]{wtc}} uses to smooth
#' the wavelet transforms. Further information regarding this method of computation
#' can be accessed through the references section.
#'
#'
#'
#' @author Alvaro Chao-Ecija
#'
#' @references
#' Keissar K, Maestri R, Pinna GD, La Rovere MT, Gilad O. Non-invasive
#' baroreflex sensitivity assessment using wavelet transfer function-based
#' time-frequency analysis. Physiol Meas. 2010 ;31(7):1021-36.
#'
#' @import biwavelet
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' TransferFun <- TransferFunCWT(Data)
TransferFunCWT <- function(data,
                           HF = 0.4,
                           LF = 0.15,
                           VLF = 0.04,
                           chosen.dj = 1 / 20,
                           demean = TRUE,
                           smooth = TRUE,
                           alpha = FALSE) {
  if (!is.data.frame(data))
    data <- as.data.frame(data)
  if (demean) {
    for (n in 2:ncol(data)) {
      data[, n] <- data[, n] - mean(data[, n])
    }
  }
  time <- data[, 1]
  WTransform.x <-
    biwavelet::wt(
      cbind(time, data[, 3]),
      dj = chosen.dj,
      s0 = 1 / (HF + 0.1),
      max.scale = 1 / (VLF - 0.01)
    )
  WTransform.y <-
    biwavelet::wt(
      cbind(time, data[, 2]),
      dj = chosen.dj,
      s0 = 1 / (HF + 0.1),
      max.scale = 1 / (VLF - 0.01)
    )
  Smoothed <-
    SmoothTransforms(WTransform.x, WTransform.y, chosen.dj)
  sm.WTransform.x <- Smoothed$sm.WTransform.x
  sm.WTransform.y <- Smoothed$sm.WTransform.y
  XWTransform <- Smoothed$XWTransform
  sm.XWTransform <- Smoothed$sm.XWTransform
  Cospectrum <-
    Re(sm.XWTransform / sqrt(sm.WTransform.x * sm.WTransform.y))
  Quadrature <-
    Im(sm.XWTransform / sqrt(sm.WTransform.x * sm.WTransform.y))
  Coherence <-
    abs(sm.XWTransform) ^ 2 / (sm.WTransform.x * sm.WTransform.y)
  Phase <- atan2(Quadrature, Cospectrum)
  if (!smooth) {
    sm.WTransform.x <- abs(WTransform.x$wave) ^ 2
    sm.WTransform.y <- abs(WTransform.y$wave) ^ 2
    sm.XWTransform <- XWTransform
  }
  if (!alpha) {
    TransferFun <- sm.XWTransform / sm.WTransform.x
  } else {
    TransferFun <- sqrt(sm.WTransform.y / sm.WTransform.x)
  }
  return(
    list(
      TransferFun = TransferFun,
      Coherence = Coherence,
      Freqs = 1 / WTransform.x$period,
      Cone = WTransform.x$coi,
      Time = data[, 1],
      HF = HF,
      LF = LF,
      VLF = VLF,
      type = "brs_cwt",
      Cospectrum = Cospectrum,
      Quadrature = Quadrature,
      Phase = Phase,
      Scales = WTransform.x$scale
    )
  )
}



#########################################################
# AUXILIARY FUNCTIONS                                   #
#########################################################




# Private function to smooth the wavelet transforms for the computation of the BRS.
# This function is based on the routine that the biwavelet function wtc uses to
# calculate the wavelet coherence, using the biwavelet function smooth.wavelet.
# This function adapts this routine for the computation of the BRS (for more
# information, please check the references section at the top of this document)
SmoothTransforms <- function(x, y, chosen.dj = 1 / 20) {
  N <- nrow(x$wave)
  M <- ncol(x$wave)
  inverse_scales <- matrix(rep(1 / t(x$scale), M), ncol = M,
                           nrow = N)
  XWTransform <- (Conj(x$wave) * y$wave)
  sm.WTransform.x <- biwavelet::smooth.wavelet(inverse_scales *
                                                 (abs(x$wave) ^ 2), x$dt, chosen.dj, x$scale)
  sm.WTransform.y <- biwavelet::smooth.wavelet(inverse_scales *
                                                 (abs(y$wave) ^ 2), x$dt, chosen.dj, x$scale)
  sm.XWTransform <- biwavelet::smooth.wavelet(inverse_scales *
                                                XWTransform, x$dt, chosen.dj, x$scale)
  return(
    list(
      sm.WTransform.x = sm.WTransform.x,
      sm.WTransform.y = sm.WTransform.y,
      sm.XWTransform = sm.XWTransform,
      XWTransform = XWTransform
    )
  )
}