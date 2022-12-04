


InterpolateData <- function(x, f = 4){
                   IntFunRR <- splinefun(x$Time, x$RR, 
                       method = "monoH.FC", ties = "ordered")
                   IntFunSBP <- splinefun(x$Time, x$SBP, 
                       method = "monoH.FC", ties = "ordered")
                  # if(!is.null(x$DBP)) IntFunDBP <- splinefun(x$Time, x$DBP, 
                  #    method = "monoH.FC", ties = "ordered")
                  # IntFunPP <- splinefun(x$Time, x$PP, 
                  #     method = "monoH.FC", ties = "ordered")
                   Time = seq(x$Time[1], x$Time[NROW(x$Time)], 1/f)
                   return(list(Time = Time, RR = IntFunRR(Time), 
                          SBP = IntFunSBP(Time))) #, DBP = IntFunDBP(Time),
                             #PP = IntFunPP(Time)))
}