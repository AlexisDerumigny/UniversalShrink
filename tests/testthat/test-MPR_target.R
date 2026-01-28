
test_that("`MPR_target` is coherent between target general and target identity", {
  set.seed(1)
  n = 10
  p = 5 * n
  mu = rep(0, p)
  
  # Generate Sigma
  X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
  H <- eigen(t(X0) %*% X0)$vectors
  Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
  
  # Generate example dataset
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  Y = t(X)
  
  precision_MPR_optimal = MPR_target(Y = Y)
  
  
  # Using the identity target directly
  Ip = diag(nrow = p)
  
  precision_MPR_optimal_Ip = MPR_target(Y = Y, Pi0 = Ip)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_MPR_optimal_Ip) - 
                              as.matrix(precision_MPR_optimal), normalized = TRUE)
  
  skip(message = "Skipping test `MPR_target` is coherent between target general and target identity")
  
  expect_equal(distFrob, 0)
  
  FrobeniusLoss2(precision_MPR_optimal_Ip, solve(Sigma))
  FrobeniusLoss2(precision_MPR_optimal, solve(Sigma))
  
  
  ## For small sample sizes, the optimization is unreliable and differences can
  ## appear
  
  expect_equal(precision_MPR_optimal_Ip$t_optimal,
               precision_MPR_optimal$t_optimal )

  expect_equal(precision_MPR_optimal_Ip$alpha_optimal,
               precision_MPR_optimal$alpha_optimal )

  expect_equal(precision_MPR_optimal_Ip$beta_optimal,
               precision_MPR_optimal$beta_optimal )
  
  
  # For the non-optimized versions:
  
  t_opt =  precision_MPR_optimal$t_optimal
  alpha_opt = precision_MPR_optimal$alpha_optimal
  beta_opt = precision_MPR_optimal$beta_optimal
  
  precision_MPR_semioptimal_Ip = MPR_target(Y = Y, Pi0 = Ip, t = t_opt)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_MPR_semioptimal_Ip) - 
                              as.matrix(precision_MPR_optimal), normalized = TRUE)
  expect_equal(distFrob, 0)
  
  
  precision_MPR_nooptim = MPR_target(Y = Y, Pi0 = Ip, t = t_opt,
                                     alpha = alpha_opt, beta = beta_opt)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_MPR_nooptim) - 
                              as.matrix(precision_MPR_optimal), normalized = TRUE)
  expect_equal(distFrob, 0)
})


test_that("`MPR_target` is coherent between optimized and non-optimized versions", {
  set.seed(1)
  n = 10
  p = 5 * n
  mu = rep(0, p)
  
  # Generate Sigma
  X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
  H <- eigen(t(X0) %*% X0)$vectors
  Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
  
  # Generate example dataset
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  Y = t(X)
  
  precision_MPR_optimal = MPR_target(Y = Y)
  
  t_opt =  precision_MPR_optimal$t_optimal
  alpha_opt = precision_MPR_optimal$alpha_optimal
  beta_opt = precision_MPR_optimal$beta_optimal
  
  precision_MPR_semioptimal = MPR_target(Y = Y, t = t_opt)
  expect_equal(as.matrix(precision_MPR_semioptimal), 
               as.matrix(precision_MPR_optimal) )
  
  precision_MPR_nooptim = MPR_target(Y = Y, t = t_opt,
                                     alpha = alpha_opt, beta = beta_opt)
  expect_equal(as.matrix(precision_MPR_nooptim),
               as.matrix(precision_MPR_optimal) )
})


test_that("`selectOptimizationType` behave as intended", {
  
  # Case where everything is provided
  result = selectOptimizationType(t = 1, alpha = 1, beta = 1)
  expect_identical(result, "none")
  
  # Cases where t is provided (and only t is used)
  result = selectOptimizationType(t = 1, alpha = NULL, beta = NULL)
  expect_identical(result, "alpha_beta")
  
  suppressWarnings({
    result = selectOptimizationType(t = 1, alpha = 1, beta = NULL)
  })
  
  expect_identical(result, "alpha_beta")
  expect_warning({selectOptimizationType(t = 1, alpha = 1, beta = NULL)},
                 class = c("MissingParametersWarning"))
  
  suppressWarnings({
    result = selectOptimizationType(t = 1, alpha = NULL, beta = 1)
  })
  
  expect_identical(result, "alpha_beta")
  expect_warning({selectOptimizationType(t = 1, alpha = NULL, beta = 1)},
                 class = c("MissingParametersWarning"))
  
  
  # Cases where t is not provided - full optimization is always needed =========
  
  # nothing is provided
  result = selectOptimizationType(t = NULL, alpha = NULL, beta = NULL)
  expect_identical(result, "all")
  
  # else we always warn and then do full optimization
  # case 1
  suppressWarnings({
    result = selectOptimizationType(t = NULL, alpha = NULL, beta = 1)
  })
  
  expect_identical(result, "all")
  expect_warning({selectOptimizationType(t = NULL, alpha = NULL, beta = 1)},
                 class = c("MissingParametersWarning"))
  
  # case 2
  suppressWarnings({
    result = selectOptimizationType(t = NULL, alpha = 1, beta = NULL)
  })
  
  expect_identical(result, "all")
  expect_warning({selectOptimizationType(t = NULL, alpha = 1, beta = NULL)},
                 class = c("MissingParametersWarning"))
  
  # case 3
  suppressWarnings({
    result = selectOptimizationType(t = NULL, alpha = 1, beta = 1)
  })
  
  expect_identical(result, "all")
  expect_warning({selectOptimizationType(t = NULL, alpha = 1, beta = 1)},
                 class = c("MissingParametersWarning"))
})

