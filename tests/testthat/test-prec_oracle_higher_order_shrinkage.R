test_that("m = 0 oracle uses the true covariance moments", {
  set.seed(4)
  
  X = matrix(rnorm(60), nrow = 12, ncol = 5)
  Sigma = diag(c(1, 2, 3, 4, 5))
  
  q1 = tr(Sigma) / ncol(X)
  q2 = tr(Sigma %*% Sigma) / ncol(X)
  expected = (q1 / q2) * diag(ncol(X))
  
  for (name in c("Moore-Penrose", "ridge", "MPR")) {
    result = prec_oracle_higher_order_shrinkage(
      X = X,
      m = 0,
      Sigma = Sigma,
      nameEstimator = name
    )
    
    expect_equal(as.matrix(result), expected)
    expect_equal(as.numeric(coef(result)), q1 / q2)
  }
})

