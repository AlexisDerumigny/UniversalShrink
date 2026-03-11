
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
#' @param method_inversion a character string of length 1 describing the
#' numerical inversion method to be used. Possible choices are \itemize{
#'   \item \code{"solve"}: use the \code{solve} function;
#'   
#'   \item \code{"Woodbury"}: TODO
#'   
#'   \item \code{"auto"}: this is the default. It chooses \code{"solve"} when
#'   the concentration ratio is smaller than 1 and else \code{"Woodbury"}.
#' }
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
#' \doi{10.1214/25-AOS2602}. ArXiv: \doi{10.48550/arXiv.2403.15792}.
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
#'   precision_ridge = ridge(X, t = t)
#'   
#'   cat("t = t, loss =", LossFrobenius2(precision_ridge, Sigma = Sigma), "\n")
#' }
#' 
#' 
#' @export
ridge <- function (X, centeredCov = TRUE, t, verbose = 0,
                   method_inversion = "auto")
{
  call_ = match.call()
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (method_inversion == "auto"){
    if (c_n < 1){
      method_ = "solve"
    } else {
      method_ = "Woodbury"
    }
  } else {
    method_ = method_inversion
  }
  
  if (method_ == "solve"){
    # Sample covariance matrix
    S <- cov_with_centering(X = X, centeredCov = centeredCov)
    
    iS_ridge <- solve(S + t * Ip)
  } else if (method_ == "Woodbury"){
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
    centralTerm = t(X_adjusted) %*% solve(XtX_over_n + t * In_adj) %*% X_adjusted / n_adjusted
    iS_ridge = (Ip / t) - centralTerm / t
  }
  
  result = list(
    estimated_precision_matrix = iS_ridge,
    t = t,
    n = n,
    p = p,
    centeredCov = centeredCov,
    method = "Ridge",
    method_ridge_inversion = method_,
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}



