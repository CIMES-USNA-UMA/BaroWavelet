




#' Simulation of cardiovascular data
#'
#' Data simulation with different baroreflex behaviors. Based on the algorithm described by
#' de Boer and Karemaker (please consult the "references" section for more information)
#'
#' @param use.noise Boolean. Add gaussian noise to the simulation
#' @param seed Constant for reproducibility purposes
#' @param v Vagal (instantaneous) component of the simulation. Default is 9
#' @param s Vector of sympathetic (lagged) components of the simulation
#'
#' @return A data frame containing timestamps and simulated RR and SBP recordings
#'
#' @author Alvaro Chao-Ecija
#'
#' @references
#' de Boer RW, Karemaker JM. Cross-Wavelet Time-Frequency Analysis Reveals Sympathetic Contribution to Baroreflex Sensitivity as Cause of Variable Phase Delay Between Blood Pressure and Heart Rate. Front Neurosci. 2019 Jul 9;13:694. Available from: https://www.frontiersin.org/articles/10.3389/fnins.2019.00694/full
#'
#' @export
#'
#' @examples
#' Data <- DataSimulation()
DataSimulation <-
  function(use.noise = TRUE,
           seed = 1,
           v = 9,
           s = c(0.25, 1.5, 2.75, 2, 0.75))
  {
    set.seed(seed)
    mRR <- rnorm(2500, mean = 1000, sd = 5)
    time <- cumsum(mRR / 1000) - mRR[1] / 1000
    dt <- c(1, diff(time))
    SBP <- 5 * sin(2 * pi * 0.1 * time) + 5 * sin(2 * pi *
                                                    0.25 * time) + use.noise * rnorm(2500, mean = 0,
                                                                                     sd = 2)
    RR <-
      c(rep(v, 500), rep(0, 500), rep(v, 500), rep(0, 500), rep(v, 500)) * SBP
    for (n in c(6:1000, 2006:2500))
      RR[n] <- RR[n] + sum(SBP[(n - 1):(n - 5)] * s)
    set.seed(seed + 1)
    RR <- RR + use.noise * rnorm(2500, mean = 0, sd = 5)
    SBP <- SBP + 120
    mRR <- rep(1000, 2500)
    RR <- RR + mRR
    return(data.frame(Time = time, RR = RR, SBP = SBP))
  }
