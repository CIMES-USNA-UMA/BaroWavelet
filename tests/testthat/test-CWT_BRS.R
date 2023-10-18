

Data <- InterpolateData(DataSimulation(), f = 1)


test_that("Transfer function should be computed without error", {
  expect_true(is.list(TransferFunCWT(Data)))
})

test_that("Alpha index should be computed without error", {
  expect_true(is.list(TransferFunCWT(Data, alpha = TRUE)))
})


test_that("Non smoothing should give a warning", {
  expect_warning(is.list(TransferFunCWT(Data, smooth = FALSE)))
})
