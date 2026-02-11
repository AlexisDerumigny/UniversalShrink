
#' Plain Ridge (Tikhonov) estimator without shrinkage
#' 
#' This function computes a 'plain' Ridge (Tikhonov) estimator of the precision matrix
#' \eqn{\mathbf{S}_n^-(t)=(\mathbf{S}_n+t\mathbf{I}_p)^{-1}},
#' where \eqn{\mathbf{S}_n} is the sample covariance matrix, \eqn{\mathbf{I}_p} is 
#' \eqn{p}-dimensional identity matrix and \eqn{t>0} is a given penalty parameter.
#' 
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param t parameter of the estimation.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the precision matrix, of class
#' `EstimatedPrecisionMatrix`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
#' 
#' 
#' @examples
#' 
#' n = 100
#' p = 5 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' 
#' for (t in c(0.2, 0.5, 1)){
#'   precision_ridge = ridge_no_shrinkage(X, t = t)
#' 
#'   cat("t = t, loss =", FrobeniusLoss2(precision_ridge, Sigma = Sigma), "\n")
#' }
#' 
#' 
#' @export
ridge_no_shrinkage <- function (X, centeredCov = TRUE, t, verbose = 0){
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  iS_ridge <- solve(S + t * Ip)
  
  result = list(
    estimated_precision_matrix = iS_ridge,
    t = t
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

