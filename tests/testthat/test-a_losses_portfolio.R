test_that("LossFrobenius2 also works for nmeric vectors", {
  Sigma = diag(1:5)
  
  X <- MASS::mvrnorm(n = 3, mu = rep(0,5), Sigma = Sigma)
  weights3 = GMV_Moore_Penrose(X)
  
  trueWeights = rowSums(solve(Sigma)) / sum(solve(Sigma))
  
  l1 = LossFrobenius2(weights3, trueWeights, normalized = FALSE)
  l2 = LossFrobenius2(trueWeights, weights3, normalized = FALSE)
  l3 = LossFrobenius2(as.numeric(weights3), trueWeights, normalized = FALSE)
  l4 = LossFrobenius2(trueWeights, as.numeric(weights3), normalized = FALSE)
  
  expect_all_equal(c(l2, l3, l4), l1)
})
