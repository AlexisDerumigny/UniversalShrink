test_that("as.matrix errors for invalid x type", {
  x = list()
  class(x) <- "PrecisionMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
  
  class(x) <- "CovarianceMatrix"
  expect_error(as.matrix(x), class = "UniversalShrinkError")
})


test_that("get_coefficient extracts a scalar alpha", {
  object <- list(alpha = 0.4)
  
  result <- get_coefficient(object)
  
  expect_identical(result, c(alpha = 0.4))
})

test_that("get_coefficient extracts higher-order alpha", {
  object <- list(alpha = matrix(c(0.1, 0.2, 0.3), ncol = 1))
  
  result <- get_coefficient(object)
  
  expect_identical(
    result,
    c(alpha_0 = 0.1, alpha_1 = 0.2, alpha_2 = 0.3)
  )
})

test_that("get_coefficient extracts alpha and beta", {
  object <- list(alpha = 0.4, beta = 0.6)
  
  result <- get_coefficient(object)
  
  expect_identical(
    result,
    c(alpha = 0.4, beta = 0.6)
  )
})

test_that("get_coefficient extracts optimal alpha and beta", {
  object <- list(
    alpha_optimal = 0.25,
    beta_optimal = 0.75
  )
  
  result <- get_coefficient(object)
  
  expect_identical(
    result,
    c(alpha = 0.25, beta = 0.75)
  )
})

test_that("ordinary and optimal fields can be resolved independently", {
  object <- list(
    alpha_optimal = 0.25,
    beta = 0.75
  )
  
  result <- get_coefficient(object)
  
  expect_identical(result, c(alpha = 0.25, beta = 0.75))
})

test_that(
  "get_coefficient returns an empty vector when no coefficient is stored", {
  object <- list(method = "Estimator without coefficients")
  
  expect_identical(get_coefficient(object), numeric(0))
})

test_that("get_coefficient rejects beta without alpha", {
  object <- list(beta = 0.5)
  
  expect_error(
    get_coefficient(object),
    class = "InvalidCoefficientRepresentationError"
  )
})


test_that("get_coefficient rejects vector alpha combined with beta", {
  object <- list(
    alpha = c(0.1, 0.2, 0.3),
    beta = 0.4
  )
  
  expect_error(
    get_coefficient(object),
    class = "InvalidCoefficientRepresentationError"
  )
})