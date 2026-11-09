

#' Finite sample nonlinear oracle estimator of the precision matrix
#' 
#' This function computes the finite-sample rotation-equivariant
#' oracle estimator of the precision matrix given by
#' \deqn{
#' \hat{\boldsymbol{\Pi}}_{\mathrm{NLoracle}}
#' =\sum_{i=1}^p d_i^{\mathrm{or}}\mathbf{u}_i\mathbf{u}_i^\top,
#' }
#' where the oracle eigenvalues are defined as
#' \deqn{
#' d_i^{\mathrm{or}}
#' =\frac{\mathbf{u}_i^\top\boldsymbol{\Sigma}\mathbf{u}_i}
#' {\mathbf{u}_i^\top\boldsymbol{\Sigma}^2\mathbf{u}_i},
#' \quad i=1,\ldots,p.
#' }
#' Here, \eqn{\mathbf{u}_i} denotes the \eqn{i}-th eigenvector
#' of the sample covariance matrix \eqn{\mathbf{S}_n}, and
#' \eqn{\boldsymbol{\Sigma}} is the population covariance matrix.
#' This
#'  procedure ensures that the loss 
#'  \eqn{||\hat{\boldsymbol{\Pi}}_{\mathrm{NLoracle}}\boldsymbol{\Sigma}-\mathbf{I}_p||^2_F}
#' is minimized with probability one for a finite sample. The estimator requires
#' knowledge of the population covariance
#' matrix and is therefore used as an oracle benchmark.
#'
#' @template param-X
#' @param Sigma true covariance matrix
#' 
#' @inheritParams cov_with_centering
#' 
#' @template returns-EstimatedPrecisionMatrix
#' 
#' @examples
#' 
#' p = 10
#' Sigma = diag(1:p)
#' n = 5
#' 
#' X <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
#' precision_oracle = inverse_NL_oracle(X, Sigma = Sigma, centeredCov = TRUE)
#' 
#' @export
inverse_NL_oracle <- function(X, Sigma, centeredCov = TRUE, verbose = 0)
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

