
#' Simulation of cardiovascular data
#'
#' Data simulation with different baroreflex behaviors. Based on the algorithm described by
#' de Boer and Karemaker (please consult the "references" section for more information)
#'
#' @param use.noise Boolean. Add gaussian noise to the simulation
#' @param seed Constant for reproducibility purposes
#' @param varies Select if the simulation will have variable brs or variable noise distributions
#' @param v Vagal (instantaneous) component of the simulation. Default is 9
#' @param s Vector of sympathetic (lagged) components of the simulation
#' @param sds1 Vector of standard deviations for IBI noise sources. If the noise distributions remain
#'             constant, only the first value of the vector will be considered
#' @param sds2 Vector of standard deviations for SBP noise sources. If the noise distributions remain
#'             constant, only the first value of the vector will be considered
#' @param N Length of intervals of interest, only used when noise distributions vary 
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
           varies = c("brs", "noise"),
           v = 9,
           s = c(0.25, 1.5, 2.75, 2, 0.75),
           sds1 = c(2, 4, 10, 20, 2, 4, 10, 20),
           sds2 = c(5, 10, 50, 70, 5, 10, 50, 70),
           N = 500)
  {
    GetSimulatedEV <-
      function(f,
               v,
               s) {
        # For future feature
        lags <- 0:(NROW(c(v,s)) - 1)
        angle <- -lags * f * 360 
        phase <- angle * pi / 180
        x <- y <- c(v,s)
        for (n in 1:NROW(c(v,s))) {
          x[n] <- cos(phase[n]) * c(v,s)[n]
          y[n] <- sin(phase[n]) * c(v,s)[n]
        }
        brs <- c(sum(x), sum(y))
        mbrs <- sqrt(brs[1] ^ 2 + brs[2] ^ 2)
        angle <- atan2(brs[2], brs[1]) * 180 / pi
        symp <- c(sum(x[-1]), sum(y[-1]))
        msymp <- sqrt(symp[1] ^ 2 + symp[2] ^ 2)
        symp_a <- atan2(symp[2], symp[1]) * 180 / pi
        return(
          list(
            x  = x,
            y = y,
            angles = phase * 180 / pi,
            symp = symp,
            msymp = msymp,
            asymp = symp_a,
            brs = brs,
            mbrs = mbrs,
            angle = angle
          )
        )
      }
    mRR <- mSBP <- mRR2 <-  NULL
    set.seed(seed)
    varies <- match.arg(varies)
    if (varies == "brs") {
      mRR <- rnorm(2500, mean = 1000, sd = 5)
    } else if (varies == "noise") {
      mRR <- rnorm(N * NROW(sds1), mean = 1000, sd = 5)
    }
    time <- cumsum(mRR / 1000) - mRR[1] / 1000
    if (NROW(v) > 1)
      stop("Use a single number for v")
    if (varies == "brs") {
      SBP <- 5 * sin(2 * pi * 0.1 * time) + 5 * sin(2 * pi *
                                                      0.25 * time) + use.noise * rnorm(2500, mean = 0,
                                                                                       sd = sds1[1])
      RR <-
        c(rep(v, 500), rep(0, 500), rep(v, 500), rep(0, 500), rep(v, 500)) * SBP
      for (n in c(6:1000, 2006:2500))
        RR[n] <- RR[n] + sum(SBP[(n - 1):(n - 5)] * s)
      set.seed(seed + 1)
      RR <- RR + use.noise * rnorm(2500, mean = 0, sd = sds2[1])
      SBP <- SBP + 120
      RR <- RR + 1000
    } else if (varies == "noise") {
      if(NROW(sds1) != NROW(sds2)) stop("Vectors for noise sd values must have equal length.")
      for (i in 1:NROW(sds1)) {
        set.seed(seed + 1)
        mRR2 <- c(mRR2, rnorm(N, mean = 0, sd = sds2[i]))
        set.seed(seed + 2)
        mSBP <- c(mSBP, rnorm(N, mean = 0, sd = sds1[i]))
      }
      SBP <- 5 * sin(2 * pi * 0.1 * time) + 5 * sin(2 * pi *
                                                      0.25 * time) + mSBP
      RR <- rep(v, NROW(mSBP)) * SBP
      for (n in c(6:NROW(mSBP)))
        RR[n] <- RR[n] + sum(SBP[(n - 1):(n - 5)] * s)
      SBP <- SBP + 120
      RR <- RR + mRR2 + 1000
    }
    return(data.frame(Time = time, RR = RR, SBP = SBP))
  }



