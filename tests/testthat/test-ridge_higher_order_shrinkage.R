
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


