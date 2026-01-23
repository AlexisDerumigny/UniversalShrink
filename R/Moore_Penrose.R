

#' Moore-Penrose inverse of the sample covariance matrix
#' 
#' Having a centered (by default) or non-centered observation matrix 
#' \eqn{\tilde{\mathbf{Y}}_n} with the corresponding sample covariance matrix 
#' \eqn{\mathbf{S}_n}. The Moore-Penrose inverse of this sample
#' covariance matrix, denoted by \eqn{\mathbf{S}^+_n} is computed by the following
#'  formula
#' \deqn{
#' \mathbf{S}_n^+=\left(\frac{1}{n}\tilde{\mathbf{Y}}_n\tilde{\mathbf{Y}}_n^\top\right)^+
#' =\frac{1}{\sqrt{n}}\tilde{\mathbf{Y}}_n\left(\frac{1}{n}
#' \tilde{\mathbf{Y}}_n^\top\tilde{\mathbf{Y}}_n\right)^{-2}
#' \frac{1}{\sqrt{n}}\tilde{\mathbf{Y}}_n^\top.
#' } See, the beginning of the proof of Theorem 2.1 in Bodnar and Parolya (2026)
#' for the details how the centering was done.
#' 
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param centeredCov Boolean: do we center (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix, i.e. the inverse of the covariance matrix).
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
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
#' # Convert to matrix class for computations
#' iS_MP = as.matrix(iS_MP)
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
    
    n_adjusted = n - 1
  } else {
    S <- Y %*% t(Y)/n
    Ytilde <- Y
    
    n_adjusted = n
  }
  # Inverse companion covariance
  iYtilde <- solve(t(Ytilde) %*% Ytilde / n_adjusted)
  
  # Moore-Penrose inverse
  if (n_adjusted < p){
    iS_MP <- MASS::ginv(S)
  } else if (n_adjusted > p){
    # Explicit expression which could outperform `MASS::ginv`
    iS_MP <- Ytilde %*% iYtilde %*% iYtilde %*% t(Ytilde) / (n-1)
  } else {
    stop("This estimator is not defined for p = n - 1 in the centered case,",
         "and for p = n in the non-centered case.")
  }
  
  result = list(
    estimated_precision_matrix = iS_MP
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


