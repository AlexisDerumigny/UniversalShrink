
# Compute the matrix M for the higher-order shrinkage of the Moore-Penrose-Ridge
# estimator
# 
# @return a list with
# - a square matrix M of size (m + 1)
# - a vector hm of size (m + 1)
# - the estimator of v
compute_M_t_MPR <- function(m, c_n, S_t_inverse, q1, q2, t, verbose)
{
  s_and_v = compute_sv_ridge(m = m + 1, # we need an additional power compared
                                        # to the ridge
                             c_n = c_n,
                             S_t_inverse = S_t_inverse, q1 = q1, q2 = q2, t = t,
                             verbose = verbose)
  s_ridge = s_and_v$s
  v = s_and_v$v
  
  # Computation of the adapted version of s  ===================================
  s = matrix(nrow = 2 * m, ncol = 2)
  for (j in 1:(2*m)){
    s[j, ] = 0
    # FIXME
  }
  
  # Computation of M  ==========================================================
  s2 = c(q2, s[, 2])
  
  M <- s2[1:(m+1)]
  for (j in 2:(m+1))
  {
    M <- cbind(M, s2[j:(m+j)])
  }
  
  # Computation of hm  =========================================================
  hm = c(q1, s[1:m, 1])
  
  if (verbose > 0){
    cat("Estimation of hm:\n")
    print(hm)
    cat("\n")
  }
  
  return(list(M = M, hm = hm, v = v))
}



#' Moore-Penrose-Ridge higher order shrinkage
#' 
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' @param t penalization parameter
#' 
#' @inheritParams cov_with_centering
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
#' t = 10
#' 
#' precision_MPR_Cent = MPR_shrinkage(X, centeredCov = TRUE, t = t)
#' precision_MPR_NoCent = MPR_shrinkage(X, centeredCov = FALSE, t = t)
#'
#' LossFrobenius2(precision_MPR_Cent, Sigma = Sigma)
#' LossFrobenius2(precision_MPR_NoCent, Sigma = Sigma)
#' 
#' for (m in 1:3){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       MPR_higher_order_shrinkage(X, m = m, centeredCov = TRUE, t = t)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       MPR_higher_order_shrinkage(X, m = m, centeredCov = FALSE, t = t)
#'       
#'   print(LossFrobenius2(precision_higher_order_shrinkage_Cent, Sigma = Sigma))
#'   
#'   print(LossFrobenius2(precision_higher_order_shrinkage_NoCent, Sigma = Sigma))
#' }
#' 
#' 
#' precision_higher_order_shrinkage_Cent = MPR_higher_order_shrinkage(X, m = 1, t = 100)
#' 
#' precision_shrinkage_identity_semioptimal_Cent = MPR_shrinkage(X, t = 100)
#'       
#' precision_shrinkage_general_semioptimal_Cent = MPR_shrinkage(X, t = 100, Pi0 = diag(p))
#'
#' precision_higher_order_shrinkage_Cent$alpha
#' precision_shrinkage_identity_semioptimal_Cent$beta_optimal
#' precision_shrinkage_identity_semioptimal_Cent$alpha_optimal
#' 
#' precision_shrinkage_general_semioptimal_Cent$beta_optimal
#' precision_shrinkage_general_semioptimal_Cent$alpha_optimal
#' 
#' precision_higher_order_shrinkage_Cent$M
#' precision_shrinkage_identity_semioptimal_Cent$M
#' 
#' precision_higher_order_shrinkage_Cent$hm
#' precision_shrinkage_identity_semioptimal_Cent$hm
#' 
#' LossFrobenius2(precision_higher_order_shrinkage_Cent, Sigma = Sigma)
#' LossFrobenius2(precision_shrinkage_identity_semioptimal_Cent, Sigma = Sigma)
#' 
#' 
#' @export
#' 
MPR_higher_order_shrinkage <- function(X, m, centeredCov = TRUE, t, verbose = 0)
{
  if (verbose > 0){
    cat("Starting `MPR_higher_order_shrinkage`...\n")
  }
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  if (verbose > 0){
    cat("*  t = ", t, "\n")
    cat("*  m = ", m, "\n")
  }
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  iS_ridge <- solve(S + t * Ip)
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  if (verbose > 0){
    cat("Starting values: \n")
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n\n")
  }
  
  estimatedM = compute_M_t_MPR(m = m, c_n = c_n, q1 = q1, q2 = q2,
                               S_t_inverse = iS_ridge,
                               t = t, verbose = verbose)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  if (verbose > 0){
    cat("Optimal alpha: \n")
    print(alpha)
  }
  
  estimated_precision_matrix = alpha[1] * Ip
  power_MPR = Ip
  
  for (k in 1:m){
    power_MPR = power_MPR %*% MPR_estimator
    estimated_precision_matrix = estimated_precision_matrix + 
      alpha[k + 1] * power_MPR
  }
  
  result = list(
    estimated_precision_matrix = estimated_precision_matrix,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

