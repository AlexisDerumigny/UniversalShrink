test_that("NormFrobenius2 is coherent between implementations", {
  
  NormFrobenius2_transpose <- function(M, normalized){
    FrobNorm2 = tr( M %*% t(M) )
    if (normalized){
      p = ncol(M)
      return (FrobNorm2 / p)
    } else {
      return (FrobNorm2)
    }
  }
  
  M = matrix(1:9, nrow = 3, ncol = 3)
  
  for (normalized in c(TRUE, FALSE))
  {
    old = NormFrobenius2_transpose(M = M, normalized = normalized)
    new = NormFrobenius2(M = M, normalized = normalized)
    expect_equal(new, old)
  }
})


test_that("NormFrobenius2 handles vectors", {
  x <- c(1, 2, 3)
  
  expect_equal(
    NormFrobenius2(x, normalized = FALSE),
    sum(x^2)
  )
  
  expect_equal(
    NormFrobenius2(x, normalized = TRUE),
    mean(x^2)
  )
})


test_that("vector Frobenius loss is symmetric", {
  x <- c(1, 2, 3)
  y <- c(3, 2, 1)
  
  expect_equal(
    LossFrobenius2(x, y, normalized = FALSE),
    LossFrobenius2(y, x, normalized = FALSE)
  )
})


