test_that("as.matrix errors for invalid x type", {
  x = list()
  class(x) <- "PrecisionMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
  
  class(x) <- "CovarianceMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
})
