test_that("Both method_invM = 'solve' and 'recursive' give the same result for small m", {
  
  n = 10
  
  for (p in c(n / 2, 2 * n)) {
    for (m in 1:2){
      # Remark: we start to observe differences for m = 5
      
      p = 2 * n
      mu = rep(0, p)
      
      # Generate Sigma
      X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
      H <- eigen(t(X0) %*% X0)$vectors
      Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
      
      # Generate example dataset
      X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
      
      precision_solve = Moore_Penrose_higher_order_shrinkage(
        X, m = m, centeredCov = TRUE, method_invM = "solve")
      
      precision_recursive = Moore_Penrose_higher_order_shrinkage(
        X, m = m, centeredCov = TRUE, method_invM = "recursive")
      
      expect_equal(precision_solve$invM, precision_recursive$invM)
  
    }
  }
})
