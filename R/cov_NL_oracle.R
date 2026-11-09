

#' Finite sample nonlinear oracle estimator of the covariance matrix
#'
#' This function computes the finite-sample rotation-equivariant
#' oracle estimator of the covariance matrix given by
#' \deqn{
#' \hat{\boldsymbol{\Sigma}}_{\mathrm{NLoracle}}
#' =\sum_{i=1}^p d_i^{\mathrm{or}}\mathbf{u}_i\mathbf{u}_i^\top,
#' }
#' where the oracle eigenvalues are defined as
#' \deqn{
#' d_i^{\mathrm{or}}
#' =\mathbf{u}_i^\top\boldsymbol{\Sigma}\mathbf{u}_i,
#' \quad i=1,\ldots,p.
#' }
#' Here, \eqn{\mathbf{u}_i} denotes the \eqn{i}-th eigenvector
#' of the sample covariance matrix \eqn{\mathbf{S}_n}, and
#' \eqn{\boldsymbol{\Sigma}} is the population covariance matrix.
#' The estimator minimizes the Frobenius loss within the class
#' of rotation-equivariant estimators for the realized sample
#' eigenvectors. This procedure ensures that the loss 
#' \eqn{||\hat{\boldsymbol{\Sigma}}_{\mathrm{NLoracle}}-\boldsymbol{\Sigma}||^2_F}
#' is minimized with probability one for a finite sample.
#'
#' @references
#' Ledoit and Peche (2011).
#' Eigenvectors of some large sample covariance matrix ensembles.
#' Probability Theory and Related Fields, 151, 233-264.
#'
#' Ledoit and Wolf (2012).
#' Nonlinear shrinkage estimation of large-dimensional
#' covariance matrices.
#' The Annals of Statistics, 40(2), 1024-1060.
#'
#' @template param-X
#' @param Sigma true covariance matrix
#' 
#' @inheritParams cov_with_centering
#' 
#' @template returns-EstimatedCovarianceMatrix
#' 
#' @examples
#' p = 10
#' Sigma = diag(1:p)
#' n = 5
#' 
#' # True mean is zero
#' 
#' X <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = Sigma)
#' cov_usual = cov_with_centering(X, centeredCov = TRUE)
#' cov_oracle = cov_NL_oracle(X, Sigma = Sigma, centeredCov = TRUE)
#' 
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

