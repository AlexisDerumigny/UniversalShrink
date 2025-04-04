

#' Quadratic inverse shrinkage
#' 
#' 
#' This function estimates the covariance matrix as a given data set using
#' quadratic inverse shrinkage.
#' 
#' @param Y data matrix (rows are observations, columns are features).
#' 
#' @returns the estimator of the covariance matrix
#' (a `p` by `p` matrix).
#' 
#' @references 
#' Ledoit, O., & Wolf, M. (2022).
#' Quadratic shrinkage for large covariance matrices.
#' Bernoulli, 28(3), 1519-1547.
#' \link{https://doi.org/10.3150/20-BEJ1315}
#' 
#' @examples
#' p = 200
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
#' estimatedCov_sample = cov(X)
#' estimatedCov_shrink = quadratic_inverse_shrinkage(X)
#' 
#' # We now compare the distance between the true and both estimators.
#' mean((eigen(Sigma)$values - eigen(estimatedCov_sample)$values)^2)
#' mean((eigen(Sigma)$values - eigen(estimatedCov_shrink)$values)^2)
#' 
#' @export
quadratic_inverse_shrinkage <- function(Y, k = NA) {
  # Y: raw data matrix
  # k: degree of freedom adjustment (default: 1)
  
  N <- nrow(Y) # Sample size
  p <- ncol(Y) # Matrix dimension
  
  # Default settings
  if (is.na(k) || is.null(k)) {
    Y <- scale(Y, center = TRUE, scale = FALSE) # Demean the raw data matrix
    k <- 1
  }
  
  n <- N - k  # Adjust effective sample size
  c <- p / n  # Concentration ratio
  
  sample <- (t(Y) %*% Y) / n  # Sample covariance matrix
  
  eig_decomp <- eigen(sample, symmetric = TRUE) # Spectral decomposition
  lambda <- sort(eig_decomp$values)  # Sorted eigenvalues (ascending)
  u <- eig_decomp$vectors[, order(eig_decomp$values)] # Corresponding eigenvectors
  
  # COMPUTE Quadratic-Inverse Shrinkage estimator of the covariance matrix
  h <- min(c^2, 1/c^2)^0.35 / p^0.35  # Smoothing parameter
  
  invlambda <- 1 / lambda[max(1, p - n + 1):p] # Inverse of non-null eigenvalues
  
  Lj <- matrix(rep(invlambda, each = min(p, n)), ncol = min(p, n)) # Like 1/lambda_j
  Lj_i <- Lj - t(Lj) # (1/lambda_j) - (1/lambda_i)
  
  theta <- rowMeans(Lj * Lj_i / (Lj_i^2 + h^2 * Lj^2)) # Smoothed Stein shrinker
  Htheta <- rowMeans(Lj * (h * Lj) / (Lj_i^2 + h^2 * Lj^2)) # Conjugate term
  Atheta2 <- theta^2 + Htheta^2 # Squared amplitude
  
  if (p <= n) {
    # Case where sample covariance matrix is not singular
    delta <- 1 / ((1 - c)^2 * invlambda + 2 * c * (1 - c) * invlambda * theta +
                    c^2 * invlambda * Atheta2) # Optimally shrunk eigenvalues
  } else {
    # Case where sample covariance matrix is singular
    delta0 <- 1 / ((c - 1) * mean(invlambda)) # Shrinkage of null eigenvalues
    delta <- c(rep(delta0, p - n), 1 / (invlambda * Atheta2))
  }
  
  deltaQIS <- delta * (sum(lambda) / sum(delta)) # Preserve trace
  
  sigmahat <- u %*% diag(deltaQIS) %*% t(u) # Reconstruct covariance matrix
  
  return(sigmahat)
}

