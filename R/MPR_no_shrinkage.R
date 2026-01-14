
#' Moore-Penrose-Ridge
#' 
#' This function computes
#' \deqn{\widehat{\Sigma^{-1}}^{ridge}_t = (S + t I_p)^{-1} - t * (S + t I_p)^{-2}},
#' where \eqn{S} is the sample covariance matrix and \eqn{t} is a given parameter.
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param t parameter of the estimation.
#' 
#' @returns the estimator of the precision matrix, of class
#' `EstimatedPrecisionMatrix`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2024).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \link{https://doi.org/10.48550/arXiv.2403.15792}
#' 
#' 
#' @examples
#' 
#' n = 10
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
#' Y = t(X)
#' Ip = diag(nrow = p)
#' for (t in c(0.01, 0.2, 0.5, 1, 2)){
#'   precision_MPR_Cent = 
#'     MPR_no_shrinkage(Y = t(X), centeredCov = TRUE, t = t)
#'   
#'   Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
#'   S = Y %*% Jn %*% t(Y) / (n-1)
#'   precision_ridge = solve(S + t * Ip)
#'   precision_MPR_Cent_alternative_expression = precision_ridge %*% S %*% precision_ridge
#'   precision_MP = Moore_Penrose(Y = t(X), centeredCov = TRUE)
#' 
#'   cat("t = ", t,", Moore-Penrose  , loss =", 
#'     FrobeniusLoss2(precision_MP, Sigma = Sigma, type = "precision"), "\n")
#'   cat("t = ", t,", ridge          , loss =", 
#'     FrobeniusLoss2(precision_ridge, Sigma = Sigma, type = "precision"), "\n")
#'   cat("t = ", t,", MPR            , loss =", 
#'     FrobeniusLoss2(precision_MPR_Cent, Sigma = Sigma, type = "precision"), "\n")
#'   cat("t = ", t,", MPR alternative, loss =", 
#'     FrobeniusLoss2(precision_MPR_Cent_alternative_expression, Sigma = Sigma,
#'                    type = "precision"), "\n")
#' }
#' 
#' 
#' @export
MPR_no_shrinkage <- function (Y, centeredCov, t){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
  } else {
    S <- Y %*% t(Y)/n
  }
  
  
  iS_ridge <- solve(S + t * Ip)
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  result = list(
    estimated_precision_matrix = MPR_estimator,
    t = t
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

