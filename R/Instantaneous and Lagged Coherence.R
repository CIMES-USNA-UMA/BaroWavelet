






CalculateCausalCouplings <- function(tf){
  Cospectrum <- tf$Cospectrum
  Quadrature <- tf$Quadrature
  Instantaneous <- abs(Cospectrum)^2
  Lagged <- (abs(Quadrature)^2) / (1 - Instantaneous)
  output <- list(Instantaneous = Instantaneous, Lagged = Lagged)
  return(output)
}
