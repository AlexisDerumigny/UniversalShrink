
test_that("non-finite Bell polynomials produce a structured error", {
  Bell_polynomials <- matrix(
    c(1, 2, NA, 4, Inf, 6),
    nrow = 2,
    ncol = 3
  )
  
  condition <- tryCatch(
    {
      stop_non_finite_bell_polynomials(
        Bell_polynomials = Bell_polynomials,
        m = 2,
        c_n = 1.5,
        q1 = 0.1,
        q2 = 0.2,
        v = c(1, 2, 3)
      )
    },
    error = function(cnd) cnd
  )
  
  # Check that an error was actually raised and caught.
  expect_s3_class(condition, "error")
  
  # Check the package-specific error classes.
  expect_s3_class(condition, "BellPolynomialNumericalError")
  expect_s3_class(condition, "NumericalError")
  expect_s3_class(condition, "UniversalShrinkError")
  
  # Check that the complete matrix is stored in the condition.
  expect_identical(condition$Bell_polynomials, Bell_polynomials)
  
  # Check the number of non-finite indices
  expect_identical(condition$n_non_finite, 2L)
})


test_that(
  "ridge higher-order shrinkage with m = 1 equals linear ridge shrinkage",
{
  set.seed(123)
  
  n <- 12
  p <- 4
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  t_value <- 0.8
  identity <- diag(p)
  
  ridge_matrix <- as.matrix( ridge(X, centeredCov = TRUE, t = t_value) )
  
  higher_order <- ridge_higher_order_shrinkage(
    X, m = 1, centeredCov = TRUE, t = t_value)
  
  higher_order_coefficients <- coef(higher_order)
  
  expected_higher_order <-
    higher_order_coefficients["alpha_0"] * identity +
    higher_order_coefficients["alpha_1"] * ridge_matrix
  
  expect_s3_class(higher_order, "EstimatedPrecisionMatrix")
  expect_named(higher_order_coefficients, c("alpha_0", "alpha_1") )
  
  expect_equal(as.matrix(higher_order), expected_higher_order,
               tolerance = 1e-10)
  
  linear_shrinkage <- ridge_shrinkage(X, centeredCov = TRUE,
                                      t = t_value)
  
  linear_shrinkage_coefficients <- coef(linear_shrinkage)
  
  expected_linear_shrinkage <-
    linear_shrinkage_coefficients["alpha"] * ridge_matrix +
    linear_shrinkage_coefficients["beta"] * identity
  
  expect_s3_class(linear_shrinkage, "EstimatedPrecisionMatrix")
  expect_named(linear_shrinkage_coefficients, c("alpha", "beta") )
  
  expect_equal(as.matrix(linear_shrinkage), expected_linear_shrinkage,
               tolerance = 1e-10)
  
  expect_equal(unname(higher_order_coefficients["alpha_0"]),
                      unname(linear_shrinkage_coefficients["beta"]),
               tolerance = 1e-10)
  
  expect_equal(unname(higher_order_coefficients["alpha_1"]),
                      unname(linear_shrinkage_coefficients["alpha"]),
               tolerance = 1e-10)
  
  expect_equal(as.matrix(higher_order), as.matrix(linear_shrinkage),
               tolerance = 1e-10)
} )
