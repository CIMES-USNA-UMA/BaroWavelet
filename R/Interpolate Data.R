
#' Interpolate data
#'
#' Interpolates data up to a certain sample frequency
#' @param x Multivariate time series to be interpolated
#' @param f Sample frequency. Default is 4 Hz

#' @return An interpolated time series.
#'
#' @author Alvaro Chao-Ecija
#'
#'
#' @export
#'
#' @examples
#' Data <- InterpolateData(DataSimulation(), f = 1)
#' 
InterpolateData <- function(x, f = 4){
                   IntFunRR <- splinefun(x$Time, x$RR, 
                       method = "monoH.FC", ties = "ordered")
                   IntFunSBP <- splinefun(x$Time, x$SBP, 
                       method = "monoH.FC", ties = "ordered")
                  # Code for future potential implementation:
                  # if(!is.null(x$DBP)) IntFunDBP <- splinefun(x$Time, x$DBP, 
                  #    method = "monoH.FC", ties = "ordered")
                  # IntFunPP <- splinefun(x$Time, x$PP, 
                  #     method = "monoH.FC", ties = "ordered")
                   Time = seq(x$Time[1], x$Time[NROW(x$Time)], 1/f)
                   return(list(Time = Time, RR = IntFunRR(Time), 
                          SBP = IntFunSBP(Time))) 
                   # Code for future potential implementation:
                          #, DBP = IntFunDBP(Time),
                             #PP = IntFunPP(Time)))
}