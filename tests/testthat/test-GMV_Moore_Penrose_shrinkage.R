
test_that("GMV_Moore_Penrose_shrinkage is coherent for eq weighted targets", {
  set.seed(1)
  n = 50
  p = 2 * n
  mu = rep(0, p)
  
  # Generate Sigma
  X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
  H <- eigen(t(X0) %*% X0)$vectors
  Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
  
  # Generate example dataset
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  GMV_MP_shrinkage_Cent_eq = 
    GMV_Moore_Penrose_shrinkage_eq(X = X, centeredCov = TRUE, verbose = 0)
    
  GMV_MP_shrinkage_Cent_gen = 
    GMV_Moore_Penrose_shrinkage_general(X = X, centeredCov = TRUE, verbose = 0)
  
  expect_equal(as.numeric(GMV_MP_shrinkage_Cent_gen),
               as.numeric(GMV_MP_shrinkage_Cent_eq))
})


test_that("default and explicit equal-weight targets agree", {
  X0 <- matrix(c(1, 0,
                 0, 1,
                 1, 1,
                 2, -1), ncol = 2, byrow = TRUE)

  for (X in list(X0, t(X0))) {
    p <- ncol(X)

    for (centered in c(TRUE, FALSE)) {
      default_result <- as.numeric(
        GMV_Moore_Penrose_shrinkage(X, centeredCov = centered))
      explicit_result <- as.numeric(
        GMV_Moore_Penrose_shrinkage(
          X, centeredCov = centered, b = rep(1 / p, p)))

      expect_equal(explicit_result, default_result)
    }
  }
})
