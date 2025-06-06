

#' Moore-Penrose inverse of the sample covariance matrix
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param centeredCov Boolean: do we center (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @examples
#' 
#' n = 50
#' p = 2 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
#' 
#' iS_MP = Moore_Penrose(Y = t(X), centeredCov = TRUE)
#' 
#' # Sample covariance matrix
#' Y = t(X)
#' Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
#' S <- Y %*% Jn %*% t(Y) / (n-1)
#' 
#' # We test iS_MP %*% S %*% iS_MP == iS_MP
#' sum((iS_MP %*% S %*% iS_MP - iS_MP)^2)
#' 
#' # We test S %*% iS_MP %*% S == S
#' sum((S %*% iS_MP %*% S - S)^2)
#' 
#' # We test S %*% iS_MP == t(S %*% iS_MP)
#' sum((S %*% iS_MP - t(S %*% iS_MP))^2)
#' 
#' # We test iS_MP %*% S == t(iS_MP %*% S)
#' sum((iS_MP %*% S - t(iS_MP %*% S))^2)
#' 
#' @export
Moore_Penrose <- function(Y, centeredCov)
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    # We remove the last eigenvector because the eigenvalues are sorted
    # in decreasing order.
    Hn = eigen(Jn)$vectors[, -n]
    Ytilde = Y %*% Hn
    
    # Inverse companion covariance
    iYtilde <- solve(t(Ytilde) %*% Ytilde / (n-1) )
    
    # Moore-Penrose inverse
    iS_MP <- Ytilde %*% iYtilde %*% iYtilde %*% t(Ytilde) / (n-1)
    
  } else {
    S <- Y %*% t(Y)/n
    
    # Inverse companion covariance
    iY <- solve(t(Y) %*% Y / n)
    
    # Moore-Penrose inverse
    iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
  }
  
  return (iS_MP)
}


