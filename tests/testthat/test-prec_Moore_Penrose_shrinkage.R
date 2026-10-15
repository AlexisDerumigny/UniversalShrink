test_that("`Moore_Penrose_shrinkage` uses a supplied target", {
  X <- matrix(c(1, 0,
                0, 1,
                1, 1,
                2, -1), ncol = 2, byrow = TRUE)
  p <- ncol(X)
  Pi0 <- diag(c(1, 2))

  result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE, Pi0 = Pi0))
  
  identity_result <- as.matrix(
    Moore_Penrose_shrinkage(X, centeredCov = FALSE, Pi0 = diag(p)))

  expect_false(isTRUE(all.equal(result, identity_result)))
})

test_that(paste0("`Moore_Penrose_shrinkage_general_plarge` and\n",
                 "`Moore_Penrose_shrinkage_identity_plarge` give coherent results"), {
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
     Moore_Penrose_shrinkage_general_plarge(X = X, centeredCov = TRUE)
     
  precision_MoorePenrose_NoCent = 
     Moore_Penrose_shrinkage_general_plarge(X = X, centeredCov = FALSE)
  
  precision_MoorePenrose_Cent_id =
    Moore_Penrose_shrinkage_identity_plarge(X = X, centeredCov = TRUE)
  
  precision_MoorePenrose_NoCent_id = 
    Moore_Penrose_shrinkage_identity_plarge(X = X, centeredCov = FALSE)
  
  expect_equal(precision_MoorePenrose_Cent, precision_MoorePenrose_Cent_id)
  expect_equal(precision_MoorePenrose_NoCent, precision_MoorePenrose_NoCent_id)
})


test_that("default and explicit identity targets give the same result", {
  X <- matrix(c(1, 0,
                0, 1,
                1, 1,
                2, -1), ncol = 2, byrow = TRUE)
  p <- ncol(X)

  for (centered in c(TRUE, FALSE)) {
    default_result <- as.matrix(
      Moore_Penrose_shrinkage(X, centeredCov = centered))
    
    identity_result <- as.matrix(
      Moore_Penrose_shrinkage(X, centeredCov = centered, Pi0 = diag(p)))
    
    expect_equal(identity_result, default_result)
  }
  
  for (centered in c(TRUE, FALSE)) {
    default_result <- as.matrix(
      Moore_Penrose_shrinkage(t(X), centeredCov = centered))
    
    identity_result <- as.matrix(
      Moore_Penrose_shrinkage(t(X), centeredCov = centered, Pi0 = diag(4)))
    
    expect_equal(identity_result, default_result)
  }
})


test_that("Moore_Penrose_shrinkage errors at the non-centered boundary p = n", {
  n <- 5L
  p <- n
  X <- diag(n)
  
  condition <- expect_error(Moore_Penrose_shrinkage(X, centeredCov = FALSE),
                            class = "UndefinedMoorePenroseError")
  
  expect_identical(condition$p, p)
  expect_identical(condition$n, n)
  expect_identical(condition$centeredCov, FALSE)
  expect_identical(condition$estimatorName, "Moore_Penrose_shrinkage")
  
  expect_match(conditionMessage(condition), "p = n", fixed = TRUE)
})


test_that("Moore_Penrose_shrinkage errors at the centered boundary p = n - 1", {
  n <- 6L
  p <- n - 1L
  X <- diag(n)[, seq_len(p), drop = FALSE]
  
  condition <- expect_error(Moore_Penrose_shrinkage(X, centeredCov = TRUE),
                            class = "UndefinedMoorePenroseError"
  )
  
  expect_identical(condition$p, p)
  expect_identical(condition$n, n)
  expect_identical(condition$centeredCov, TRUE)
  expect_identical(condition$estimatorName, "Moore_Penrose_shrinkage")
  
  expect_match(conditionMessage(condition), "p = n - 1", fixed = TRUE)
})


test_that("Moore_Penrose_shrinkage accepts p = n for centered covariance", {
  n <- 5L
  p <- n
  X <- diag(n)
  
  expect_no_error(result <- Moore_Penrose_shrinkage(X, centeredCov = TRUE))
  
  estimate <- as.matrix(result)
  
  expect_s3_class(result, "EstimatedPrecisionMatrix")
  expect_identical(dim(estimate), c(p, p))
  expect_true(all(is.finite(estimate)))
  expect_equal(estimate, t(estimate), tolerance = 1e-12)
})

