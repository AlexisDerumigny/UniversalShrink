test_that("LossFrobenius2 also works for numeric vectors", {
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


test_that("LossFrobenius2 computes the expected loss for numeric vectors", {
  x <- c(0.2, 0.3, 0.5)
  y <- c(0.1, 0.4, 0.5)
  
  result <- LossFrobenius2(
    x,
    otherPortfolioWeights = y,
    normalized = FALSE
  )
  
  expected <- sum((y - x)^2)
  
  expect_equal(result, expected)
})


test_that("LossFrobenius2 computes the normalized loss", {
  x <- c(0.2, 0.3, 0.5)
  y <- c(0.1, 0.4, 0.5)
  
  result <- LossFrobenius2(
    x,
    otherPortfolioWeights = y,
    normalized = TRUE
  )
  
  expected <- mean((y - x)^2)
  
  expect_equal(result, expected)
})


test_that("LossFrobenius2 rejects vectors with different lengths", {
  x <- c(0.2, 0.3, 0.5)
  y <- c(0.4, 0.6)
  
  expect_error(
    LossFrobenius2(x, otherPortfolioWeights = y),
    class = "IncompatibleVectorLengthsError"
  )
})


test_that("LossFrobenius2 rejects non-numeric portfolio weights", {
  x <- c(0.2, 0.3, 0.5)
  y <- c("0.1", "aaa", "0.5")
  
  suppressWarnings({
    expect_error(
      LossFrobenius2(x, otherPortfolioWeights = y),
      class = "InvalidArgumentError"
    )
  })
})


test_that("LossFrobenius2 rejects non-finite portfolio weights", {
  x <- c(0.2, 0.3, 0.5)
  y <- c(0.1, NA_real_, 0.9)
  
  expect_error(
    LossFrobenius2(x, otherPortfolioWeights = y),
    class = "InvalidArgumentError"
  )
})

