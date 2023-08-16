

test_that("Lack of data should return error", {
  expect_error(InterpolateData())
})

test_that("Wrong frequency value should return an error", {
  Data <- DataSimulation()
  expect_error(InterpolateData(Data, f = "f"))
})

test_that("Wrong data set length should return an error", {
  Data <- data.frame(Time = 0, RR = 950, SBP = 120)
  expect_error(InterpolateData(Data))
})

test_that("Interpolate a signal with dt equal to f should be fine", {
  set.seed(1)
  Data <-
    data.frame(
      Time = 0:100,
      RR = rnorm(101, 1000),
      SBP = rnorm(101, 120)
    )
  Int_Data <- InterpolateData(Data, f = 1)
  expect_equal(Int_Data$RR, Data$RR)
})



