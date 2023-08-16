


Data <- InterpolateData(DataSimulation(), f = 1)


test_that("Function should work without error", {
  expect_true(is.list(AlphaIndexDWT(Data, wv = "d8", error = 0.01)))
})

test_that("HRV data should be computed when required", {
  expect_true(is.list(AlphaIndexDWT(
    Data,
    wv = "d8",
    error = 0.01,
    hrv = TRUE
  )$HRV))
})


test_that("Function should give error when no data available", {
  expect_error(AlphaIndexDWT("Data", wv = "d8", error = 0.01))
})

test_that("Function should give error when no wavelet available", {
  expect_error(AlphaIndexDWT(Data, wv = "wv", error = 0.01))
})

test_that("Function should give error when no error available", {
  expect_error(AlphaIndexDWT(Data, wv = "d8", error = "error"))
})



