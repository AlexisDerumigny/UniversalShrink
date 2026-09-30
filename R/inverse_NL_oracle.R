

#' Oracle estimator of the precision matrix
#'
#' @param X data matrix (rows are observations, columns are features).
#' @param Sigma true covariance matrix
#' 
#' @inheritParams cov_with_centering
#' 
#' @export
inverse_NL_oracle <- function(X, Sigma, centeredCov = TRUE, verbose = 0)
{
  call_ = match.call()
  
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
  
  Sigma2 = Sigma %*% Sigma
  
  new_eigenvalues = diag(t(U) %*% Sigma %*% U) / diag(t(U) %*% Sigma2 %*% U)
  
  inverse_oracle <- U %*% diag(new_eigenvalues) %*% t(U)
  
  
  result = new_estimated_precision_matrix(
    estimate = inverse_oracle,
    n = n,
    p = p,
    centeredCov = centeredCov,
    method = "inverse_NL_oracle",
    call = call_
  )
  
  return (result)
}

