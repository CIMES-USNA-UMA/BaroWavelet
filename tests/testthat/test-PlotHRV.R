pdf(NULL)


Data <- InterpolateData(DataSimulation(), f = 1)
AlphaDWT <- AlphaIndexDWT(Data, hrv = TRUE)



test_that("HRV should be plotted as LF/HF ratio", {
  expect_null(PlotHRV(
    AlphaDWT$HRV,
    AlphaDWT$Time,
    newPlot = FALSE,
    use.ggplot = FALSE
  ))
})


test_that("HRV should be plotted showing both LF and HF", {
  expect_null(
    PlotHRV(
      AlphaDWT$HRV,
      AlphaDWT$Time,
      newPlot = FALSE,
      use.ggplot = FALSE,
      ratio = FALSE
    )
  )
})


test_that("Only LF should be shown", {
  expect_null(
    PlotHRV(
      AlphaDWT$HRV,
      AlphaDWT$Time,
      newPlot = FALSE,
      use.ggplot = FALSE,
      ratio = FALSE,
      plotHF = FALSE
    )
  )
})


test_that("Using ggplot style should be fine", {
  expect_true(is.list(
    PlotHRV(
      AlphaDWT$HRV,
      AlphaDWT$Time,
      newPlot = FALSE,
      use.ggplot = TRUE,
      plotHF = FALSE
    )
  ))
})


test_that("Not using HRV should return an error", {
  expect_error(is.list(PlotHRV(NULL, AlphaDWT$Time, newPlot = FALSE)))
})

test_that("Not using a time support should return an error", {
  expect_error(is.list(PlotHRV(AlphaDWT$HRV, NULL, newPlot = FALSE)))
})
