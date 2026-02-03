
#' Ridge with target set to the identity
#' 
#' This function computes
#' \eqn{\widehat{\Sigma^{-1}}^{ridge}_t = (S + t I_p)^{-1}},
#' where \eqn{S} is the sample covariance matrix and \eqn{t} is a given parameter.
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param t parameter of the estimation.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the precision matrix, of class
#' `EstimatedPrecisionMatrix`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2024).
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
#'   precision_ridge_Cent = 
#'     ridge_no_shrinkage(Y = t(X), centeredCov = TRUE, t = t)
#' 
#'   cat("t = t, loss =", FrobeniusLoss2(precision_ridge_Cent, Sigma = Sigma), "\n")
#' }
#' 
#' 
#' @export
ridge_no_shrinkage <- function (Y, centeredCov = TRUE, t, verbose = 0){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  iS_ridge <- solve(S + t * Ip)
  
  result = list(
    estimated_precision_matrix = iS_ridge,
    t = t
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

