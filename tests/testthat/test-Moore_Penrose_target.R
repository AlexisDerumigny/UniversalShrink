test_that("`Moore_Penrose_target_general` and `Moore_Penrose_target_identity` give coherent results", {
  set.seed(1)
  n = 50
  p = 5 * n
  mu = rep(0, p)
  
  # Generate Sigma
  X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
  H <- eigen(t(X0) %*% X0)$vectors
  Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
  
  # Generate example dataset
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
  
  precision_MoorePenrose_Cent =
     Moore_Penrose_target_general(X = X, centeredCov = TRUE)
     
  precision_MoorePenrose_NoCent = 
     Moore_Penrose_target_general(X = X, centeredCov = FALSE)
  
  precision_MoorePenrose_Cent_id =
    Moore_Penrose_target_identity(X = X, centeredCov = TRUE)
  
  precision_MoorePenrose_NoCent_id = 
    Moore_Penrose_target_identity(X = X, centeredCov = FALSE)
  
  expect_equal(precision_MoorePenrose_Cent, precision_MoorePenrose_Cent_id)
  expect_equal(precision_MoorePenrose_NoCent, precision_MoorePenrose_NoCent_id)
})
