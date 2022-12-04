





#' @export
PreprocessData <- function(data, maxSBP = 300, minSBP = 30,
                           minbpm = 25, maxbpm = 200, use.RHRV = TRUE){
  Time <- data[[1]]
  SBP <- data$SBP
  DBP <- data$DBP
  N <- NROW(Time)
  for(n in 2:N) if(Time[n]<Time[n-1]) break
  if(!is.null(n)) Time[n:N] <- Time[n:N]*60
  N <- NROW(Time)
  for(m in 1:N) if(SBP[m] > maxSBP | SBP[m] < minSBP) Time[m] <- NA
  SBP <- SBP[!is.na(Time)]
  DBP <- DBP[!is.na(Time)]
  Time <- Time[!is.na(Time)]
  if(use.RHRV){
    RHRV_object <- PrepareTachoForRHRV(Time)
    RHRV_object <- RHRV::BuildNIHR(RHRV_object)
    RHRV_object <- RHRV::FilterNIHR(RHRV_object)
    RR <- RHRV_object$Beat$RR
    SBP <- SBP[Time %in% RHRV_object$Beat$Time]
    DBP <- DBP[Time %in% RHRV_object$Beat$Time]
    Time <- RHRV_object$Beat$Time
  } else {
    N <- NROW(Time)
    RR <- double(N)
    RR[2:N] <- diff(Time * 1000)
    RR[1] <- RR[2]
    HR <- 60000 / RR
    for(m in 1:N) if(HR[m] > maxbpm | HR[m] < minbpm) Time[m] <- NA
    RR <- RR[!is.na(Time)]
    SBP <- SBP[!is.na(Time)]
    DBP <- DBP[!is.na(Time)]
    Time <- Time[!is.na(Time)]
  }
  return(list(Time = Time, RR = RR, SBP = SBP, DBP = DBP, PP = SBP - DBP, 
              MAP = (2*SBP - DBP)/3))
}
  
  




#' @export
PrepareTachoForRHRV <- function(time){
  HRV <- list()
  HRV$Ext <- "hrv"
  HRV$FreqAnalysis <- list()
  HRV$Verbose <- FALSE
  HRV$Beat <- list()
  HRV$Beat$Time <- time
  return(HRV)
}
