
# Tests for invalid matrix estimates ===========================================
  
test_that("matrix constructors reject non-matrix estimates", {
  expect_error(
    new_estimated_precision_matrix(c(1, 2, 3), n = 2, p = 3,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_covariance_matrix(data.frame(x = 1:3), n = 2, p = 3,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
})


test_that("matrix constructors reject non-numeric matrices", {
  character_matrix <- matrix(letters[1:4], nrow = 2)
  
  expect_error(
    new_estimated_precision_matrix(character_matrix, n = 2, p = 2,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_covariance_matrix(character_matrix, n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
})


test_that("matrix constructors reject empty matrices", {
  empty_matrix <- matrix(numeric(), nrow = 0, ncol = 0)
  
  expect_error(
    new_estimated_precision_matrix(empty_matrix, n = 2, p = 0,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X))),
    "must not be empty",
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_covariance_matrix(empty_matrix, n = 2, p = 0,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must not be empty",
    class = "InvalidEstimatorObjectError"
  )
})


test_that("matrix constructors reject non-square matrices", {
  non_square_matrix <- matrix(1:6, nrow = 2, ncol = 3)
  
  expect_error(
    new_estimated_precision_matrix(non_square_matrix, n = 2, p = 2,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X))),
    "must be a square matrix",
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_covariance_matrix(non_square_matrix, n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must be a square matrix",
    class = "InvalidEstimatorObjectError"
  )
})

test_that("portfolio constructor rejects non-numeric objects", {
  expect_error(
    new_estimated_portfolio_weights(c("a", "b"), n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must be a numeric vector",
    class = "InvalidEstimatorObjectError"
  )
})


# Tests for invalid portfolio estimates  =======================================

test_that("portfolio constructor rejects matrices and arrays", {
  expect_error(
    new_estimated_portfolio_weights(matrix(c(0.4, 0.6), ncol = 1), n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must be a numeric vector",
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_portfolio_weights(array(1:8, dim = c(2, 2, 2)), n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must be a numeric vector",
    class = "InvalidEstimatorObjectError"
  )
})


test_that("portfolio constructor rejects empty vectors", {
  expect_error(
    new_estimated_portfolio_weights(numeric(), n = 2, p = 0,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    "must not be empty",
    class = "InvalidEstimatorObjectError"
  )
})


# Tests for dimension consistency  =============================================

test_that("constructors accept a consistent supplied p", {
  expect_no_error(
    new_estimated_precision_matrix(diag(3), n = 2, p = 3,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X)))
  )
  
  expect_no_error(
    new_estimated_covariance_matrix(diag(4), n = 2, p = 4,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X)))
  )
  
  expect_no_error(
    new_estimated_portfolio_weights(rep(0.2, 5), n = 2, p = 5,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X)))
  )
})


test_that("constructors reject an inconsistent supplied p", {
  expect_error(
    new_estimated_precision_matrix(diag(3), n = 2, p = 4,
                                   centeredCov = TRUE,
                                   method = "Test estimator",
                                   call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_covariance_matrix(diag(3), n = 2, p = 2,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
  
  expect_error(
    new_estimated_portfolio_weights(c(0.3, 0.7), n = 2, p = 3,
                                    centeredCov = TRUE,
                                    method = "Test estimator",
                                    call = quote(test_estimator(X))),
    class = "InvalidEstimatorObjectError"
  )
})


test_that("constructors reject malformed supplied p", {
  invalid_values <- list(NA, Inf, numeric(), c(2, 2), "2")
  
  for (invalid_p in invalid_values) {
    expect_error(
      new_estimated_precision_matrix(diag(2), n = 2, p = invalid_p,
                                     centeredCov = TRUE,
                                     method = "Test estimator",
                                     call = quote(test_estimator(X))),
      class = "InvalidEstimatorObjectError"
    )
  }
})

