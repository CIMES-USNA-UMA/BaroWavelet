Data <- InterpolateData(DataSimulation(), f = 1)
AlphaDWT <- AlphaIndexDWT(Data)
TransferFun <- TransferFunCWT(Data)

test_that("Indices should be computed with no problems", {
  expect_true(is.matrix(IndividualIndices(AlphaDWT, method = "median")))
  expect_true(is.matrix(IndividualIndices(TransferFun, method = "median")))
  expect_true(is.matrix(
    IndividualIndices(TransferFun, method = "median", use.phase = TRUE)
  ))
})


test_that("Indices should be computed with no problems #2", {
  expect_true(is.matrix(IndividualIndices(AlphaDWT, method = "mean")))
  expect_true(is.matrix(IndividualIndices(TransferFun, method = "mean")))
  expect_true(is.matrix(
    IndividualIndices(TransferFun, use.phase = TRUE, method = "mean")
  ))
})



