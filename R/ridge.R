
#' Plain Ridge (Tikhonov) estimator without shrinkage
#' 
#' This function computes a 'plain' Ridge (Tikhonov) estimator of the precision
#' matrix
#' \eqn{\mathbf{S}_n^-(t)=(\mathbf{S}_n+t\mathbf{I}_p)^{-1}},
#' where \eqn{\mathbf{S}_n} is the sample covariance matrix, \eqn{\mathbf{I}_p}
#' is the \eqn{p}-dimensional identity matrix and \eqn{t>0} is a given penalty
#' parameter.
#' 
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param t parameter of the estimation.
#' 
#' @param method_inversion a character string of length 1 describing the
#' numerical inversion method to be used. Possible choices are \itemize{
#'   \item \code{"solve"}: compute the inverse directly using
#'   \code{\link[base]{solve}}.
#'
#'   \item \code{"woodbury"}: compute the inverse using the Woodbury matrix
#'   identity. This can be more efficient when the sample size is smaller than
#'   the matrix dimension.
#'
#'   \item \code{"auto"}: choose \code{"solve"} when the adjusted concentration
#'   ratio is smaller than one, and \code{"woodbury"} otherwise.
#' }
#' 
#' @inheritParams cov_with_centering
#' 
#' 
#' @template returns-EstimatedPrecisionMatrix
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.1214/25-AOS2602}. ArXiv: \doi{10.48550/arXiv.2403.15792}.
#' 
#' 
#' @examples
#' 
#' n = 100
#' p = 2 * n
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
#'   precision_ridge = ridge(X, t = t)
#'   
#'   cat("t = t, loss =", LossFrobenius2(precision_ridge, Sigma = Sigma), "\n")
#' }
#' 
#' 
#' @export
ridge <- function (X, centeredCov = TRUE, t, verbose = 0,
                   method_inversion = c("auto", "solve", "woodbury"))
{
  call_ = match.call()
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  method_inversion <- match.arg(method_inversion)
  if (method_inversion == "auto"){
    c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov,
                         verbose = verbose - 1)
    if (c_n < 1){
      method_ = "solve"
    } else {
      method_ = "woodbury"
    }
  } else {
    method_ = method_inversion
  }
  
  if (method_ == "solve"){
    # Sample covariance matrix
    S <- cov_with_centering(X = X, centeredCov = centeredCov)
    
    iS_ridge <- solve(S + t * Ip)
  } else if (method_ == "woodbury"){
    if (centeredCov){
      n_adjusted = n - 1
      
      Jn = diag(n) - matrix(1/n, nrow = n, ncol = n)
      eig_decomp = eigen(Jn)
      U = eig_decomp$vectors
      # we delete the last column
      Hn = U[, -n]
      X_adjusted = t(Hn) %*% X
    } else {
      n_adjusted = n
      X_adjusted = X
    }
    In_adj = diag(n_adjusted)
    
    XtX_over_n = X_adjusted %*% t(X_adjusted) / n_adjusted
    
    centralTerm = 
      t(X_adjusted) %*% solve(XtX_over_n + t * In_adj) %*% X_adjusted
    centralTerm = centralTerm / n_adjusted
    
    iS_ridge = (Ip / t) - centralTerm / t
  }
  
  result = new_estimated_precision_matrix(
    estimate = iS_ridge,
    n = n,
    p = p,
    centeredCov = centeredCov,
    t = t,
    method = "Ridge",
    method_ridge_inversion = method_,
    call = call_
  )
  
  return (result)
}



