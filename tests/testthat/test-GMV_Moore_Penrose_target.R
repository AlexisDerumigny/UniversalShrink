
test_that("GMV_Moore_Penrose_target is coherent for eq weighted targets", {
  set.seed(1)
  n = 50
  p = 2 * n
  mu = rep(0, p)
  
  # Generate Sigma
  X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
  H <- eigen(t(X0) %*% X0)$vectors
  Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
  
  # Generate example dataset
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
  
  GMV_MP_shrinkage_Cent_eq = 
    GMV_Moore_Penrose_target_eq(Y = t(X), centeredCov = TRUE, verbose = 0)
    
  GMV_MP_shrinkage_Cent_gen = 
    GMV_Moore_Penrose_target_general(Y = t(X), centeredCov = TRUE, verbose = 0)
  
  expect_equal(GMV_MP_shrinkage_Cent_gen, GMV_MP_shrinkage_Cent_eq)
})
