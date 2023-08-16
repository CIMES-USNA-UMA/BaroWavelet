
pdf(NULL)

Data <- InterpolateData(DataSimulation(), f = 1)
AlphaDWT <- AlphaIndexDWT(Data)
AlphaCWT <- TransferFunCWT(Data, alpha = TRUE)
TransferFun <- TransferFunCWT(Data)


test_that("AlphaDWT should be plotted", {
  expect_null(PlotBRS(AlphaDWT, newPlot = FALSE))
})

test_that("AlphaCWT should be plotted", {
  expect_null(PlotBRS(AlphaCWT, newPlot = FALSE))
})

test_that("TransferFun should be plotted", {
  expect_null(PlotBRS(TransferFun, newPlot = FALSE))
})

test_that("Average plot should work", {
  expect_null(PlotBRS(TransferFun, newPlot = FALSE, avg = TRUE))
})


test_that("Using ggplot style should be fine", {
  expect_true(is.list(PlotBRS(AlphaDWT, newPlot = FALSE, use.ggplot = TRUE)))
})

test_that("Wrong size axis should give error", {
  expect_error(PlotBRS(AlphaDWT, newPlot = FALSE, size.axis = "size.axis"))
})

test_that("Wrong size labels should give error", {
  expect_error(PlotBRS(AlphaDWT, newPlot = FALSE, size.axis = "size.labels"))
})

test_that("Wrong size title should give error", {
  expect_error(PlotBRS(AlphaDWT, newPlot = FALSE, size.title = "size.title"))
})

test_that("Wrong ylim should give error", {
  expect_error(PlotBRS(AlphaDWT, newPlot = FALSE, ylim = "ylim"))
})



