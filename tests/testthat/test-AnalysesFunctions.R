pdf(NULL)


test_that("A structure should be built with no problems", {
  expect_true(is.list(BuildStructure()))
})

Structure <- BuildStructure(error = 0.01)
Data <- InterpolateData(DataSimulation())

test_that("Analysis additions to Structure should work", {
  Structure <- AddAnalysis(Structure)
  expect_true(is.list(Structure$Analyses[[1]]))
})

Structure <- AddAnalysis(Structure)

test_that("Data should be added to Structure without problems", {
  Structure <-
    AddDataToAnalysis(
      Structure,
      1,
      RR = Data$RR,
      SBP = Data$SBP,
      time = Data$Time
    )
  expect_equal(Structure$Analyses[[1]]$Data[, "RR"], Data$RR)
})

Structure <-
  AddDataToAnalysis(Structure,
                    1,
                    RR = Data$RR,
                    SBP = Data$SBP,
                    time = Data$Time)

test_that("DWT analysis should work", {
  Structure <- AnalyzeBRS(Structure, 1, method = "dwt")
  expect_true(is.list(Structure$Analyses[[1]]$BRS$DWT))
})

test_that("CWT analysis should work", {
  Structure <- AnalyzeBRS(Structure, 1, method = "cwt")
  expect_true(is.list(Structure$Analyses[[1]]$BRS$CWT))
})

test_that("Both analysis methods should work", {
  Structure <- AnalyzeBRS(Structure, 1)
  expect_true(is.list(Structure$Analyses[[1]]$BRS$DWT))
  expect_true(is.list(Structure$Analyses[[1]]$BRS$CWT))
})

Structure <- AnalyzeBRS(Structure, 1)

test_that("Plotting DWT should work", {
  expect_true(is.null(
    PlotAnalyzedBRS(
      Structure,
      1,
      method = "dwt",
      use.ggplot = FALSE,
      newPlot = FALSE
    )
  ))
  
})

test_that("Plotting DWT with ggplot should work", {
  expect_true(is.list(
    PlotAnalyzedBRS(
      Structure,
      1,
      method = "dwt",
      use.ggplot = TRUE,
      newPlot = FALSE
    )
  ))
  
})

test_that("Plotting CWT should work", {
  expect_true(is.null(PlotAnalyzedBRS(
    Structure, 1, method = "cwt", newPlot = FALSE
  )))
  
})

test_that("Adding a time interval should work", {
  Structure <- AddTimeInterval(Structure)
  expect_true(is.list(Structure$IndividualIndices[[1]]))
  
})

Structure <- AddTimeInterval(Structure)

test_that("Analyzing BRS index in a given time interval should work for dwt", {
  Structure <-
    AnalyzeBRSIndices(Structure, 1, 1, c(0, 200 / 60), method = "dwt")
  expect_true(is.matrix(Structure$IndividualIndices[[1]]$HR))
  expect_true(is.numeric(Structure$IndividualIndices[[1]]$Time_DWT))
  expect_true(is.numeric(Structure$IndividualIndices[[1]]$DWT))
  expect_true(!is.numeric(Structure$IndividualIndices[[1]]$CWT))
  
})

test_that("Analyzing BRS index in a given time interval should work for cwt", {
  Structure <-
    AnalyzeBRSIndices(Structure, 1, 1, c(0, 200 / 60), method = "cwt")
  expect_true(is.numeric(Structure$IndividualIndices[[1]]$Time_CWT))
  expect_true(is.numeric(Structure$IndividualIndices[[1]]$CWT))
  
})

Structure <- AddTimeInterval(Structure)
test_that("Analyzing BRS index in a given time interval should work for both methods",
          {
            Structure <- AnalyzeBRSIndices(Structure, 1, 2, c(0, 100 / 60))
            expect_true(is.matrix(Structure$IndividualIndices[[2]]$HR))
            expect_true(is.numeric(Structure$IndividualIndices[[2]]$Time_DWT))
            expect_true(is.numeric(Structure$IndividualIndices[[2]]$Time_CWT))
            expect_true(is.numeric(Structure$IndividualIndices[[2]]$DWT))
            expect_true(is.numeric(Structure$IndividualIndices[[2]]$CWT))
            
          })

Structure <- AnalyzeBRSIndices(Structure, 1, 1, c(0, 200 / 60))
Structure <- AnalyzeBRSIndices(Structure, 1, 2, c(0, 100 / 60))

test_that("Plotting indices should work", {
  expect_true(is.null(
    PlotIndicesFromAnalysis(Structure, newPlot = FALSE, method = "dwt")
  ))
  expect_true(is.null(
    PlotIndicesFromAnalysis(Structure, newPlot = FALSE, method = "cwt")
  ))
})












