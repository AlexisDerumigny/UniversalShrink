test_that("as.matrix errors for invalid x type", {
  x = list()
  class(x) <- "EstimatedPrecisionMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
  
  class(x) <- "EstimatedCovarianceMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
})
