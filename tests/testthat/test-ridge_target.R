test_that("`ridge_target` is coherent between target general and target identity", {
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
  
  precision_ridge_optimal = ridge_target(Y = Y, verbose = 0)
  
  
  # Using the identity target directly
  Ip = diag(nrow = p)
  
  precision_ridge_optimal_Ip = ridge_target(Y = Y, Pi0 = Ip, verbose = 0)
  
  expect_equal(precision_ridge_optimal_Ip$t_optimal,
               precision_ridge_optimal$t_optimal)
  
  expect_equal(precision_ridge_optimal_Ip$alpha_optimal,
               precision_ridge_optimal$alpha_optimal)
  
  expect_equal(precision_ridge_optimal_Ip$beta_optimal,
               precision_ridge_optimal$beta_optimal)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_ridge_optimal_Ip) - 
                              as.matrix(precision_ridge_optimal), normalized = TRUE)
  
  expect_equal(distFrob, 0)
  
  FrobeniusLoss2(precision_ridge_optimal_Ip, solve(Sigma))
  FrobeniusLoss2(precision_ridge_optimal, solve(Sigma))
  
  
  ## For small sample sizes, the optimization is unreliable and differences can
  ## appear
  
  expect_equal(precision_ridge_optimal_Ip$t_optimal,
               precision_ridge_optimal$t_optimal )
  
  expect_equal(precision_ridge_optimal_Ip$alpha_optimal,
               precision_ridge_optimal$alpha_optimal )
  
  expect_equal(precision_ridge_optimal_Ip$beta_optimal,
               precision_ridge_optimal$beta_optimal )
  
  
  # For the non-optimized versions:
  
  t_opt =  precision_ridge_optimal$t_optimal
  alpha_opt = precision_ridge_optimal$alpha_optimal
  beta_opt = precision_ridge_optimal$beta_optimal
  
  precision_ridge_semioptimal_Ip = ridge_target(Y = Y, Pi0 = Ip, t = t_opt)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_ridge_semioptimal_Ip) - 
                              as.matrix(precision_ridge_optimal), normalized = TRUE)
  expect_equal(distFrob, 0)
  
  
  precision_ridge_nooptim = ridge_target(Y = Y, Pi0 = Ip, t = t_opt,
                                     alpha = alpha_opt, beta = beta_opt)
  
  distFrob = FrobeniusNorm2(as.matrix(precision_ridge_nooptim) - 
                              as.matrix(precision_ridge_optimal), normalized = TRUE)
  expect_equal(distFrob, 0)
})


test_that("`ridge_target` is coherent between optimized and non-optimized versions", {
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
  
  precision_ridge_optimal = ridge_target(Y = Y)
  
  t_opt =  precision_ridge_optimal$t_optimal
  alpha_opt = precision_ridge_optimal$alpha_optimal
  beta_opt = precision_ridge_optimal$beta_optimal
  
  precision_ridge_semioptimal = ridge_target(Y = Y, t = t_opt)
  expect_equal(as.matrix(precision_ridge_semioptimal), 
               as.matrix(precision_ridge_optimal) )
  
  precision_ridge_nooptim = ridge_target(Y = Y, t = t_opt,
                                     alpha = alpha_opt, beta = beta_opt)
  expect_equal(as.matrix(precision_ridge_nooptim),
               as.matrix(precision_ridge_optimal) )
})

