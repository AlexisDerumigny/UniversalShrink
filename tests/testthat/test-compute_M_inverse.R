test_that("check_large_m behaves as expected", {
  expect_no_warning(check_large_m(warn = TRUE,  p = 10, m = 4))
  
  expect_warning(
    check_large_m(warn = TRUE, p = 10, m = 5),
    class = "TooLarge_m_Warning"
  )
  
  expect_no_warning(check_large_m(warn = FALSE, p = 10, m = 5))
  
  expect_error(check_large_m(warn = FALSE, p = 10, m = 5.2),
               class = "InvalidArgumentError")
  expect_error(check_large_m(warn = FALSE, p = 10, m = NULL),
               class = "InvalidArgumentError")
  expect_error(check_large_m(warn = FALSE, p = 10, m = NA),
               class = "InvalidArgumentError")
})


test_that("all bona fide precision estimators agree for m = 0", {
  set.seed(1)
  
  X = matrix(rnorm(80), nrow = 10, ncol = 8)
  
  MP = Moore_Penrose_higher_order_shrinkage(
    X,
    m = 0
  )
  
  ridge = ridge_higher_order_shrinkage(
    X,
    m = 0
  )
  
  MPR = MPR_higher_order_shrinkage(
    X,
    m = 0
  )
  
  expect_equal(as.matrix(MP), as.matrix(ridge))
  expect_equal(as.matrix(MP), as.matrix(MPR))
  
  expect_equal(coef(MP), coef(ridge))
  expect_equal(coef(MP), coef(MPR))
})


test_that("t does not affect the m = 0 estimator", {
  set.seed(2)
  
  X = matrix(rnorm(80), nrow = 10, ncol = 8)
  
  ridge1 = ridge_higher_order_shrinkage(X, m = 0, t = 1)
  ridge2 = ridge_higher_order_shrinkage(X, m = 0, t = 100)
  
  MPR1 = MPR_higher_order_shrinkage(X, m = 0, t = 1)
  MPR2 = MPR_higher_order_shrinkage(X, m = 0, t = 100)
  
  expect_equal(as.matrix(ridge1), as.matrix(ridge2))
  expect_equal(as.matrix(MPR1), as.matrix(MPR2))
  expect_equal(as.matrix(ridge1), as.matrix(MPR1))
})


test_that("m = 0 covariance estimator is the scaled identity", {
  set.seed(1)
  
  X = matrix(rnorm(60), nrow = 12, ncol = 5)
  S = cov_with_centering(X, centeredCov = TRUE)
  
  result = cov_higher_order_shrinkage(
    X,
    m = 0,
    centeredCov = TRUE
  )
  
  q1 = tr(S) / ncol(X)
  expected = q1 * diag(ncol(X))
  
  expect_equal(as.matrix(result), expected)
  expect_equal(coef(result), c(alpha = q1))
  expect_equal(result$hm, q1)
  expect_equal(result$invM, 1)
})

test_that("m = 0 bona fide precision estimators are the scaled identity", {
  set.seed(2)
  
  for (centered in c(TRUE, FALSE)) {
    X = matrix(rnorm(80), nrow = 10, ncol = 8)
    
    S = cov_with_centering(X, centeredCov = centered)
    c_n = concentr_ratio(
      n = nrow(X),
      p = ncol(X),
      centeredCov = centered,
      verbose = 0
    )
    
    q1 = tr(S) / ncol(X)
    q2 = tr(S %*% S) / ncol(X) - c_n * q1^2
    expected = (q1 / q2) * diag(ncol(X))
    
    estimators = list(
      Moore_Penrose_higher_order_shrinkage(
        X,
        m = 0,
        centeredCov = centered
      ),
      ridge_higher_order_shrinkage(
        X,
        m = 0,
        centeredCov = centered
      ),
      ridge_higher_order_shrinkage(
        X,
        m = 0,
        centeredCov = centered,
        t = 2
      ),
      MPR_higher_order_shrinkage(
        X,
        m = 0,
        centeredCov = centered
      ),
      MPR_higher_order_shrinkage(
        X,
        m = 0,
        centeredCov = centered,
        t = 2
      )
    )
    
    for (result in estimators) {
      expect_equal(as.matrix(result), expected)
      expect_equal(as.numeric(coef(result)), q1 / q2)
    }
  }
})

test_that("ridge and MPR do not optimize t when m = 0", {
  set.seed(3)
  X = matrix(rnorm(80), nrow = 10, ncol = 8)
  
  controls = list(
    method = "smoothed",
    grid = c(0.1, 1, 10),
    k = 3
  )
  
  ridge_result = ridge_higher_order_shrinkage(
    X,
    m = 0,
    optimizationControls = controls
  )
  
  MPR_result = MPR_higher_order_shrinkage(
    X,
    m = 0,
    optimizationControls = controls
  )
  
  expect_null(ridge_result$result_optimization)
  expect_null(MPR_result$result_optimization)
  
  # Here t is neither optimized by the algorithm nor given by the user.
  # so the expectation is to obtain a numeric vector of length 0.
  expect_identical(get_t(ridge_result), numeric(0))
  expect_identical(get_t(MPR_result), numeric(0))
})


