
EvaluateBRS <- function(test, time_test, time_control, control = NULL,
   plot = TRUE, weight = TRUE, set.control = 1, show.CI = FALSE, plotHF = TRUE,
     plotLF = TRUE, newPlot = FALSE, title){
            if(is.null(control)){
               output <- EvaluateBRSwithoutCon(test, time_test, time_control, plot = plot, 
                 weight = weight, use.thr = TRUE, set.control = set.control, show.CI = show.CI,
                      plotHF = plotHF, plotLF = plotLF, newPlot = newPlot, title = title)
            } else {
               output<- EvaluateBRSwithCon(test, control, time_test, plot = plot, 
                 weight = weight, show.CI = show.CI)
            }
            return(output)             
}


EvaluateBRSwithCon <- function(test, control, time_test, weight = TRUE, plot = TRUE, show.CI = TRUE,
   extended  =FALSE, sig = 0.05, newPlot = TRUE){
            Estimate_test <- ExpectedValues(test, time_test, weight = weight)
            Estimate_con <- ExpectedValues(control, weight = weight)
            BRS1 <- -(Estimate_test - Estimate_con) 
            BRS2 <- Estimate_test / Estimate_con
            BRS3 <- abs(BRS1) * 100 / Estimate_con
            time_test <- time_test * 60
            select_time <- test$Time[(test$Time >= time_test[1]) &
                     (test$Time <= time_test[2])]
            select_time <- match(select_time, test$Time)
            if(weight){
                     w1 <- 1:NROW(select_time)
                     w1 <- exp(-((w1 - mean(w1))^2) / (2 * var(w1)))
                     w2 <- 1:NROW(control$Time)
                     w2 <- exp(-((w2 - mean(w2))^2) / (2 * var(w2)))
                     HF_test <- NROW(select_time) * w1 * test$HF[select_time] / 
                        sum(w1)
                     LF_test <- NROW(select_time) * w1 * test$LF[select_time] / 
                        sum(w1)
                     HF_con <- NROW(control$Time) * w2 * control$HF / sum(w2)
                     LF_con <- NROW(control$Time) * w2 * control$LF / sum(w2)
             } else {
                     HF_test <- test$HF[select_time]
                     LF_test <- test$LF[select_time]
                     HF_con <- control$HF
                     LFcon <- control$LF
             }
             Effect <- c(rep(1,NROW(control$Time)), double(
               NROW(select_time)))
             HFdata <- data.frame(HF = c(HF_con, HF_test),
               Effect = Effect)
             LFdata <- data.frame(LF = c(LF_con, LF_test),
               Effect = Effect)
             HFmodel <- lm(HF ~ 1 + Effect, data = HFdata)
             LFmodel <- lm(LF ~ 1 + Effect, data = LFdata)
             HFmodel_test <- lmtest::coeftest(HFmodel, 
               vcov. = sandwich::vcovHAC)
             LFmodel_test <- lmtest::coeftest(LFmodel, 
               vcov. = sandwich::vcovHAC)
             BRS_HF <- HFmodel$coefficients[2]
             BRS_LF <- LFmodel$coefficients[2]
             HF_con.est <- HFmodel$coefficients[1]
             LF_con.est <- LFmodel$coefficients[1]
             HFp <- HFmodel_test[2,4]
             LFp <- LFmodel_test[2,4]
             p.value <- c(HFp, LFp)
             se_HFcon <- HFmodel_test[1,2]
             se_LFcon <- LFmodel_test[1,2]
             se_BRS <- c(HFmodel_test[2,2], LFmodel_test[2,2])
             const <- abs(qt(c(sig/2, 1 - sig/2), summary(LFmodel)$df[2])[1])
             if(plot){
                PlotBRSWithCon(se_BRS[1], se_BRS[2], BRS_HF,
                   BRS_LF, const, p.value, newPlot = TRUE)
             }
             output <- round(rbind(BRS1, BRS2, BRS3, se_BRS, p.value),4)
             rownames(output)[1:4] <- c("BRS (delta)", "BRS (ratio)", "BRS (%)", 
                "SE")
             return(output)
}


PlotBRSWithCon <- function(se_HF, se_LF, BRS_HF, BRS_LF, const, pvals, newPlot = TRUE){
                if(newPlot) x11(title = "BRS")
                ci_HF = c(BRS_HF - const * se_HF, BRS_HF + const * se_HF)
                ci_LF = c(BRS_LF - const * se_LF, BRS_LF + const * se_LF)
                Max = max(c(0, ci_HF[2],ci_LF[2], BRS_LF, BRS_HF)) + 2
                Min = min(c(0, ci_HF[1], ci_LF[1], BRS_LF, BRS_HF)) - 2
                ylim = c(Min, Max)
                barplot(cbind(HF = BRS_HF, LF = BRS_LF), beside = TRUE, ylim = ylim, col = c("grey20", "grey80"), main = 
                 "Barorreflex Sensitivity")
                abline(h = 0, col = "red")
                arrows(c(1.5, 3.5), c(ci_HF[1], ci_LF[1]), c(1.5, 3.5), c(ci_HF[2], ci_LF[2]),
                   angle = 90, code = 3)
                if(pvals[1] <= 0.001){
                   text1 = "***"
                } else if(pvals[1] <= 0.01){
                   text1 = "**"
                } else if(pvals[1] <= 0.05){
                   text1 = "*"
                } else {
                   text1 = ""
                }
                if(pvals[2] <= 0.001){
                   text2 = "***"
                } else if(pvals[2] <= 0.01){
                   text2 = "**"
                } else if(pvals[2] <= 0.05){
                   text2 = "*"
                } else {
                   text2 = ""
                }
                Max1 <- max(c(0, ci_HF[2], BRS_HF)) + 1
                Max2 <- max(c(0, ci_LF[2], BRS_LF)) + 1
                text(1.5, Max1, text1, cex = 2)
                text(3.5, Max2, text2, cex = 2)
}

EvaluateBRSwithoutCon <- function(fun, time_flags1, time_flags2, plot = TRUE, 
   weight = TRUE, use.thr = TRUE, set.control = 1, show.CI = FALSE,
     plotHF = TRUE, plotLF = TRUE, newPlot = TRUE, title, sig = 0.05){
             Estimate1 <- ExpectedValues(fun, time_flags1, weight = weight)
             Estimate2 <- ExpectedValues(fun, time_flags2, weight = weight)
             BRS1 <- -(Estimate2 - Estimate1) #* 100 / Estimate1
             BRS2 <- Estimate2 / Estimate1
             BRS3 <- abs(BRS1) * 100 / Estimate1
             time_flags1 <- time_flags1 * 60
             select_time1 <- fun$Time[(fun$Time >= time_flags1[1]) &
                     (fun$Time <= time_flags1[2])]
             select_time1 <- match(select_time1, fun$Time)
             time_flags2 <- time_flags2 * 60
             select_time2 <- fun$Time[(fun$Time >= time_flags2[1]) &
                     (fun$Time <= time_flags2[2])]
             select_time2 <- match(select_time2, fun$Time)
             point <- max(select_time1)
             if(point == select_time2[1]) select_time2 <- select_time2[-1]
             interval <- c(select_time1, select_time2)
             if(weight){
                     w1 <- 1:NROW(select_time1)
                     w1 <- exp(-((w1 - mean(w1))^2) / (2 * var(w1)))
                     w2 <- 1:NROW(select_time2)
                     w2 <- exp(-((w2 - mean(w2))^2) / (2 * var(w2)))
                     HF <- c(NROW(select_time1) * w1 * fun$HF[select_time1] / 
                        sum(w1), NROW(select_time2) * w2 * fun$HF[select_time2] / 
                           sum(w2))
                     LF <- c(NROW(select_time1) * w1 * fun$LF[select_time1] / 
                        sum(w1), NROW(select_time2) * w2 * fun$LF[select_time2] / 
                           sum(w2))
             } else {
                     HF <- fun$HF[interval]
                     LF <- fun$LF[interval]
             }
             if(set.control == 1){
                Effect <- c(rep(1, NROW(select_time1)),double(NROW(select_time2)))
                Aux <- select_time1
             } else {
                Effect <- c(double(NROW(select_time1)), rep(1, NROW(select_time2)))
                Aux <- select_time2
             }
             HFdata = data.frame(HF = HF, x = fun$Time[interval], Effect = Effect)
             LFdata = data.frame(LF = LF, x = fun$Time[interval], Effect = Effect)
             HFmodel <- lm(HF ~ 1 + Effect, data = HFdata, na.action = na.omit)
             LFmodel <- lm(LF ~ 1 + Effect, data = LFdata, na.action = na.omit)
             HFtest <- car::linearHypothesis(HFmodel, c("Effect=0"), 
               vcov. = sandwich::vcovHAC)
             LFtest <- car::linearHypothesis(LFmodel, c("Effect=0"),
               vcov. = sandwich::vcovHAC)
             realHF <- lmtest::coeftest(HFmodel, vcov. = sandwich::vcovHAC)
             realLF <- lmtest::coeftest(LFmodel, vcov. = sandwich::vcovHAC)
             seHF <- realHF[1,2]
             seLF <- realLF[1,2]
             HFp <- HFtest[2,4]
             LFp <- LFtest[2,4]
             p.value <- c(HFp, LFp)
             const <- abs(qt(c(sig/2, 1 - sig/2), summary(LFmodel)$df[2])[1])
             if(plot){
                PlotBRSWithCon(se_HF, se_LF, realHF[2,1],
                   realLF[2,1], const, p.value, newPlot = TRUE)
             }    
             output <- round(rbind(BRS1, BRS2, BRS3, se_BRS, p.value),4)
             rownames(output)[1:4] <- c("BRS (delta)", "BRS (ratio)", "BRS (%)", 
                "SE")
             return(output)
}

