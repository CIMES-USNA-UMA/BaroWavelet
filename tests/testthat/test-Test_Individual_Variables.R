Data <- InterpolateData(DataSimulation(), f = 1)
AlphaDWT <- AlphaIndexDWT(Data, hrv = TRUE)
TransferFun <- TransferFunCWT(Data)

test_that("Function should work for AlphaDWT", {
  test <- TestIndBRS(AlphaDWT, c(0, 200 / 60), c(300 / 60, 450 / 60))
  expect_equal(NROW(test), 2)
})

test_that("Function should work for TransferFun", {
  test <-
    TestIndBRS(TransferFun, c(0, 200 / 60), c(300 / 60, 450 / 60))
  expect_equal(NROW(test), 2)
})

test_that("An error should be returned when equal intervals are used", {
  expect_error(TestIndBRS(TransferFun, c(0, 200 / 60), c(0 / 60, 200 / 60)))
})


test_that("Function should work for Data", {
  test <- TestIndHRandBP(Data, c(0, 200 / 60), c(300 / 60, 450 / 60))
  expect_equal(NROW(test), 2)
})


test_that("An error should be returned when equal intervals are used #2", {
  expect_error(TestIndHRandBP(Data, c(0, 200 / 60), c(0 / 60, 200 / 60)))
})


test_that("Function should work for AlphaDWT #2", {
  test <- TestIndHRV(AlphaDWT, c(0, 200 / 60), c(300 / 60, 450 / 60))
  expect_equal(NROW(test), 5)
})


test_that("An error should be returned when equal intervals are used #3", {
  expect_error(TestIndHRV(AlphaDWT, c(0, 200 / 60), c(0 / 60, 200 / 60)))
})
