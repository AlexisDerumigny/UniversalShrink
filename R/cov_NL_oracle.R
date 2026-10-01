

#' NL oracle estimator of the covariance matrix
#'
#' @param X data matrix (rows are observations, columns are features).
#' @param Sigma true covariance matrix
#' 
#' @inheritParams cov_with_centering
#' 
#' @template returns-EstimatedCovarianceMatrix
#' 
#' @export
cov_NL_oracle <- function(X, Sigma, centeredCov = TRUE, verbose = 0)
{
  call_ = match.call()
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  eigendecomposition = eigen(S)
  U = eigendecomposition$vectors
  
  NonLin_oracle <- U %*% diag(diag(t(U) %*% Sigma %*% U)) %*% t(U)
  
  result = new_estimated_covariance_matrix(
    estimate = NonLin_oracle,
    n = n,
    p = p,
    centeredCov = centeredCov,
    method = "cov_NL_oracle",
    call = call_
  )
  
  return (result)
}

