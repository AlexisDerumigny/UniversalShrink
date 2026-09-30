test_that("Moore_Penrose errors at the non-centered boundary p = n", {
  n <- 5L
  p <- n
  X <- diag(n)
  
  condition <- expect_error(Moore_Penrose(X, centeredCov = FALSE),
                            class = "UndefinedMoorePenroseError")
  
  expect_identical(condition$p, p)
  expect_identical(condition$n, n)
  expect_identical(condition$centeredCov, FALSE)
  expect_identical(condition$estimator, "Moore_Penrose")
  
  expect_match(conditionMessage(condition), "p = n", fixed = TRUE)
})


test_that("Moore_Penrose errors at the centered boundary p = n - 1", {
  n <- 6L
  p <- n - 1L
  X <- diag(n)[, seq_len(p), drop = FALSE]
  
  condition <- expect_error(Moore_Penrose(X, centeredCov = TRUE),
                            class = "UndefinedMoorePenroseError"
  )
  
  expect_identical(condition$p, p)
  expect_identical(condition$n, n)
  expect_identical(condition$centeredCov, TRUE)
  expect_identical(condition$estimator, "Moore_Penrose")
  
  expect_match(conditionMessage(condition), "p = n - 1", fixed = TRUE)
})


test_that("Moore_Penrose accepts p = n for centered covariance", {
  n <- 5L
  p <- n
  X <- diag(n)
  
  expect_no_error(result <- Moore_Penrose(X, centeredCov = TRUE))
  
  estimate <- as.matrix(result)
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  expect_identical(dim(estimate), c(p, p))
  expect_true(all(is.finite(estimate)))
  expect_equal(estimate, t(estimate), tolerance = 1e-12)
})

