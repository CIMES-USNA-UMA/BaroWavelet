
test_that("Data simulation with variable brs should work", {
  expect_true(is.list(DataSimulation()))
})

test_that("Data simulation with variable noise distributions should work", {
  expect_true(is.list(DataSimulation(varies = "noise")))
})

test_that("Wrong length for vector v should return an error", {
  expect_error(DataSimulation(v = c(9, 9)))
})

test_that("Wrong length for vector v should return an error #2", {
  expect_error(NoiseSimulation(v = c(9, 9)))
})
