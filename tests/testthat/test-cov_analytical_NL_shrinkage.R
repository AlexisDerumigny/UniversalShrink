test_that("cov_analytical_NL_shrinkage works better than sample cov", {
  set.seed(1)
  
  p = 200
  Sigma = diag(seq(1, 0.02, length.out = p))
  mu = rep(0, p)
  X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
  estimatedCov_sample = cov(X)
  estimatedCov_shrink = cov_analytical_NL_shrinkage(X)

  # We now compare the distance between the true and both estimators.
  loss_sampleF = FrobeniusLoss2(estimatedCov_sample, Sigma, type = "covariance")
  loss_shrinkF = FrobeniusLoss2(estimatedCov_shrink, Sigma)
  
  expect_lt(loss_shrinkF, loss_sampleF)

  loss_sampleE = LossEuclideanEigenvalues2(estimatedCov_sample, Sigma, type = "covariance")
  loss_shrinkE = LossEuclideanEigenvalues2(estimatedCov_shrink, Sigma)
  
  expect_lt(loss_shrinkE, loss_sampleE)
})
