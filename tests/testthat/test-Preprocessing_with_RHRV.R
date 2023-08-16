
Data <- DataSimulation()

test_that("Function should work with the data (using RHRV)", {
  expect_true(is.list(PreprocessData(Data)))
})

test_that("Function should work with the data (not using RHRV)", {
  expect_true(is.list(PreprocessData(Data, use.RHRV = FALSE)))
})
