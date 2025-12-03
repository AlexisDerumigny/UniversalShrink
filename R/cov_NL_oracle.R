

#' NL oracle estimator of the covariance matrix
#'
#' @param X data matrix (rows are observations, columns are features).
#' @param Sigma true covariance matrix
#' 
#' @inheritParams cov_with_centering
#' 
#' @export
cov_NL_oracle <- function(X, Sigma, centeredCov = TRUE, verbose = 0){
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  eigendecomposition = eigen(S)
  U = eigendecomposition$vectors
  
  NonLin_oracle <- U %*% diag(diag(t(U) %*% Sigma %*% U)) %*% t(U)
  
  result = list(
    estimated_covariance_matrix = NonLin_oracle
  )
  
  class(result) <- c("EstimatedCovarianceMatrix")
  
  return (result)
}

