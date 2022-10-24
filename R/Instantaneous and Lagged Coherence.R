






CalculateCausalCouplings <- function(tf, use.thr = TRUE, thr = 0.5){
  Cospectrum <- tf$Cospectrum
  Quadrature <- tf$Quadrature
  Instantaneous <- abs(Cospectrum)^2
  Lagged <- (abs(Quadrature)^2) / (1 - Instantaneous)
  if(use.thr){
    Coherence <- tf$Coherence
    Lagged[which(Coherence < thr)] <- NA
    Instantaneous[which(Coherence < thr)] <- NA
  }
  output <- list(Instantaneous = Instantaneous, Lagged = Lagged)
  return(output)
}
