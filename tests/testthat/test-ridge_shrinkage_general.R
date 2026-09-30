# Shared deterministic data -----------------------------------------------

set.seed(123)

n <- 12
p <- 4

X <- matrix(rnorm(n * p), nrow = n, ncol = p)

t_value <- 0.8
alpha_value <- 0.6
beta_value <- 0.3

Ip <- diag(p)

# Deliberately nonidentity target.
# This is important for detecting accidental use of Ip instead of Pi0.
Pi0 <- diag(c(1, 2, 3, 4))


# 1. Identity target, no optimization -------------------------------------

test_that("ridge_shrinkage uses supplied t, alpha, and beta with identity target", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    t = t_value,
    alpha = alpha_value,
    beta = beta_value
  )
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = t_value
  ))
  
  expected <- alpha_value * ridge_matrix + beta_value * Ip
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  expect_equal(as.matrix(result), expected, tolerance = 1e-10)
  expect_equal(result$t, t_value)
  expect_equal(result$n, n)
  expect_equal(result$p, p)
  expect_true(result$centeredCov)
})


# 2. General target, no optimization --------------------------------------

test_that("ridge_shrinkage uses supplied general target without optimization", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    Pi0 = Pi0,
    t = t_value,
    alpha = alpha_value,
    beta = beta_value
  )
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = t_value
  ))
  
  expected <- alpha_value * ridge_matrix + beta_value * Pi0
  
  expect_equal(as.matrix(result), expected, tolerance = 1e-10)
})


# 3. Identity target, fixed t and optimized alpha/beta --------------------

test_that("ridge_shrinkage optimizes alpha and beta for a fixed t", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    t = t_value
  )
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = t_value
  ))
  
  expected <- result$alpha_optimal * ridge_matrix +
    result$beta_optimal * Ip
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  
  expect_true(is.numeric(result$alpha_optimal))
  expect_length(result$alpha_optimal, 1L)
  expect_true(is.finite(result$alpha_optimal))
  
  expect_true(is.numeric(result$beta_optimal))
  expect_length(result$beta_optimal, 1L)
  expect_true(is.finite(result$beta_optimal))
  
  expect_equal(result$t, t_value)
  expect_equal(as.matrix(result), expected, tolerance = 1e-8)
})


# 4. General target, fixed t and optimized alpha/beta ---------------------

test_that("ridge_shrinkage optimizes alpha and beta for a general target", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    Pi0 = Pi0,
    t = t_value
  )
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = t_value
  ))
  
  expected <- result$alpha_optimal * ridge_matrix +
    result$beta_optimal * Pi0
  
  expect_true(is.finite(result$alpha_optimal))
  expect_true(is.finite(result$beta_optimal))
  expect_equal(as.matrix(result), expected, tolerance = 1e-8)
})


optimization_controls <- list(
  method = "smoothed",
  grid = c(0.2, 0.4, 0.8, 1.6, 3.2),
  k = 3
)


# 5. Identity target, optimize t, alpha, and beta -------------------------

test_that("ridge_shrinkage fully optimizes with identity target", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    optimizationControls = optimization_controls
  )
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  
  expect_true(is.numeric(result$t_optimal))
  expect_length(result$t_optimal, 1L)
  expect_true(is.finite(result$t_optimal))
  expect_true(result$t_optimal > 0)
  
  expect_true(is.finite(result$alpha_optimal))
  expect_true(is.finite(result$beta_optimal))
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = result$t_optimal
  ))
  
  expected <- result$alpha_optimal * ridge_matrix +
    result$beta_optimal * Ip
  
  expect_equal(as.matrix(result), expected, tolerance = 1e-8)
  expect_true(result$t_optimal %in% optimization_controls$grid)
})



# 6. General target, optimize t, alpha, and beta --------------------------

test_that("ridge_shrinkage fully optimizes with a general target", {
  result <- ridge_shrinkage(
    X,
    centeredCov = TRUE,
    Pi0 = Pi0,
    optimizationControls = optimization_controls
  )
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  
  expect_true(is.finite(result$t_optimal))
  expect_true(result$t_optimal > 0)
  expect_true(is.finite(result$alpha_optimal))
  expect_true(is.finite(result$beta_optimal))
  
  ridge_matrix <- as.matrix(ridge(
    X,
    centeredCov = TRUE,
    t = result$t_optimal
  ))
  
  expected <- result$alpha_optimal * ridge_matrix +
    result$beta_optimal * Pi0
  
  expect_equal(as.matrix(result), expected, tolerance = 1e-8)
  expect_true(result$t_optimal %in% optimization_controls$grid)
})
