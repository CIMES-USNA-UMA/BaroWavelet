

# Data simulation of physiological signals. Based on the algorithm described by
# De Boer and Karemaker in "Cross-Wavelet Time-Frequency Analysis Reveals
# Sympathetic Contribution to Baroreflex Sensitivity as Cause of Variable Phase
# Delay Between Blood Pressure and Heart Rate".


DataSimulation <- function(use.noise = TRUE){
  set.seed(1)
  SBP <- 5  *sin(2*pi*0.1*(1:2000)) + 5  *sin(2*pi*0.25*(1:2000)) +
    use.noise * rnorm(2000, mean = 0, sd = 2)
  RR <- c(rep(9,1499), rep(3, 300), rep(9, 201)) * SBP
  for(n in 7:2000) RR[n] <- RR[n] + sum(SBP[(n-2):(n-6)] * c(1.5,2.75,2,0.75,0.25))
  set.seed(2)
  RR <- RR + use.noise * rnorm(2000, mean = 0, sd = 5)
  SBP <- SBP + 120
  RR <- RR + c(rep(1000,1499), rep(700, 300), rep(1000, 201))
  return(data.frame(Time = 1:2000, RR = RR, SBP = SBP))
}
