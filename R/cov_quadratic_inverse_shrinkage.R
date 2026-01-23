

#' Quadratic inverse shrinkage (QIS)
#' 
#' 
#' This function estimates the covariance matrix as a given data set using
#' quadratic inverse shrinkage.
#' 
#' The quadratic inverse shrinkage estimator (QIS) of 
#' Ledoit and Wolf is given by 
#' \eqn{\mathbf{S}_{QIS}= 
#' \mathbf{U} \text{diag}(\hat\delta_{n,1},\ldots,\hat\delta_{n,p})\mathbf{U}^\top}
#'  with 
#' \deqn{\hat{\delta}^{-1}_{n,i}=
#' \left(1-\frac{p}{n}\right)^2 \lambda^{-1}_{n,i}+
#' 2\frac{p}{n}\left(1-\frac{p}{n}\right)
#'\lambda^{-1}_{n,i}\hat{\theta}_n(\lambda^{-1}_{n,i})+\left(\frac{p}{n}\right)^2
#'\lambda^{-1}_{n,i}\mathcal{A}^2_{\theta_n}(\lambda^{-1}_{n,i}),
#'  }
#'  for \eqn{i\in\{1,\ldots,p\}}, where 
#'  
#' \deqn{
#'\hat{\theta}_n(x)
#':=
#'  \frac{1}{p}
#'\sum_{j=1}^{p}
#'\lambda^{-1}_{n,j}
#'\frac{\lambda^{-1}_{n,j}-x}
#'{(\lambda^{-1}_{n,j}-x)^2 + h_n^2 \lambda^{-2}_{n,j}}
#'} and
#'\deqn{
#'\mathcal{A}^2_{\theta_n}(x)
#'=
#'  \left[
#'    \frac{1}{p}
#'    \sum_{j=1}^{p}
#'    \lambda^{-1}_{n,j}
#'    \frac{\lambda^{-1}_{n,j}-x}
#'    {(\lambda^{-1}_{n,j}-x)^2 + h_n^2 \lambda^{-2}_{n,j}}
#'    \right]^2
#'+
#'  \left[
#'    \frac{1}{p}
#'    \sum_{j=1}^{p}
#'    \lambda^{-1}_{n,j}
#'    \frac{h_n \lambda^{-1}_{n,j}}
#'    {(\lambda^{-1}_{n,j}-x)^2 + h_n^2 \lambda^{-2}_{n,j}}
#'    \right]^2.
#'    }
#'  The smoothing parameter \eqn{h_n} satisfies the conditions 
#'  \eqn{h_n\sim Kn^{-\alpha}} for some \eqn{K>0} and \eqn{\alpha\in(0, 2/5)}.
#'  See, Theorem 4.1 of Ledoit and Wolf (2022) for more details.
#' 
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the covariance matrix
#' (a `p` by `p` matrix).
#' TODO: update this
#' 
#' @references 
#' Ledoit, O., & Wolf, M. (2022).
#' Quadratic shrinkage for large covariance matrices.
#' Bernoulli, 28(3), 1519-1547.
#' \doi{10.3150/20-BEJ1315}
#' 
#' @examples
#' p = 200
#' n = 400
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
#' estimatedCov_sample = cov(X)
#' estimatedCov_shrink = cov_quadratic_inverse_shrinkage(X)
#' 
#' # We now compare the distance between the true and both estimators.
#' FrobeniusLoss2(estimatedCov_sample, Sigma, type = "covariance")
#' FrobeniusLoss2(estimatedCov_shrink, Sigma)
#' 
#' LossEuclideanEigenvalues2(estimatedCov_sample, Sigma, type = "covariance")
#' LossEuclideanEigenvalues2(estimatedCov_shrink, Sigma)
#' 
#' @export
cov_quadratic_inverse_shrinkage <- function(X, centeredCov = TRUE, verbose = 0) {
  
  n <- nrow(X) # Sample size
  p <- ncol(X) # Matrix dimension
  
  n_adjusted <- if (centeredCov) (n - 1) else n
  c <- concentration_ratio(n = n, p = p, centeredCov = centeredCov, 
                           verbose = verbose)
  
  # Sample covariance matrix
  sample <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  eig_decomp <- eigen(sample, symmetric = TRUE) # Spectral decomposition
  lambda <- sort(eig_decomp$values)  # Sorted eigenvalues (ascending)
  u <- eig_decomp$vectors[, order(eig_decomp$values)] # Corresponding eigenvectors
  
  # COMPUTE Quadratic-Inverse Shrinkage estimator of the covariance matrix
  h <- min(c^2, 1/c^2)^0.35 / p^0.35  # Smoothing parameter
  
  invlambda <- 1 / lambda[max(1, p - n_adjusted + 1):p] # Inverse of non-null eigenvalues
  
  # Like 1/lambda_j
  Lj <- matrix(rep(invlambda, each = min(p, n_adjusted)), ncol = min(p, n_adjusted))
  Lj_i <- Lj - t(Lj) # (1/lambda_j) - (1/lambda_i)
  
  theta <- rowMeans(Lj * Lj_i / (Lj_i^2 + h^2 * Lj^2)) # Smoothed Stein shrinker
  Htheta <- rowMeans(Lj * (h * Lj) / (Lj_i^2 + h^2 * Lj^2)) # Conjugate term
  Atheta2 <- theta^2 + Htheta^2 # Squared amplitude
  
  if (p <= n_adjusted) {
    # Theorem 4,1, Eq 4.12, page 1532-1533 of (Ledoit, O., & Wolf, M., 2022).
    # Case where sample covariance matrix is not singular
    delta <- 1 / ((1 - c)^2 * invlambda + 2 * c * (1 - c) * invlambda * theta +
                    c^2 * invlambda * Atheta2) # Optimally shrunk eigenvalues
  } else {
    # Case where sample covariance matrix is singular
    delta0 <- 1 / ((c - 1) * mean(invlambda)) # Shrinkage of null eigenvalues
    delta <- c(rep(delta0, p - n_adjusted), 1 / (invlambda * Atheta2))
  }
  
  deltaQIS <- delta * (sum(lambda) / sum(delta)) # Preserve trace
  
  sigmahat <- u %*% diag(deltaQIS) %*% t(u) # Reconstruct covariance matrix
  
  result = list(
    estimated_covariance_matrix = sigmahat
  )
  
  class(result) <- c("EstimatedCovarianceMatrix")
  
  return(result)
}

