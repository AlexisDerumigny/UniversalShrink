test_that("`Moore_Penrose_shrinkage` uses a supplied target", {
  X <- matrix(c(1, 0,
                0, 1,
                1, 1,
                2, -1), ncol = 2, byrow = TRUE)
  p <- ncol(X)
  Pi0 <- diag(p) * 2

  result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE, Pi0 = Pi0))
  
  identity_result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE, Pi0 = diag(p)))

  expect_false(isTRUE(all.equal(result, identity_result)))
})

test_that("`Moore_Penrose_shrinkage_general` and `Moore_Penrose_shrinkage_identity` give coherent results", {
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
     Moore_Penrose_shrinkage_general(X = X, centeredCov = TRUE)
     
  precision_MoorePenrose_NoCent = 
     Moore_Penrose_shrinkage_general(X = X, centeredCov = FALSE)
  
  precision_MoorePenrose_Cent_id =
    Moore_Penrose_shrinkage_identity(X = X, centeredCov = TRUE)
  
  precision_MoorePenrose_NoCent_id = 
    Moore_Penrose_shrinkage_identity(X = X, centeredCov = FALSE)
  
  expect_equal(precision_MoorePenrose_Cent, precision_MoorePenrose_Cent_id)
  expect_equal(precision_MoorePenrose_NoCent, precision_MoorePenrose_NoCent_id)
})


test_that("default and explicit identity targets give the same result in the case p < n", {
  X <- matrix(c(1, 0,
                0, 1,
                1, 1,
                2, -1), ncol = 2, byrow = TRUE)
  p <- ncol(X)

  default_result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE))
  
  identity_result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE, Pi0 = diag(p)))

  expect_equal(identity_result, default_result)
})

