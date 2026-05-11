

#' Higher order shrinkage for estimation of the covariance matrix
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns an object of class `EstimatedCovarianceMatrix` containing 
#' \itemize{
#'    \item `estimated_covariance_matrix`: the estimator of the covariance
#'    matrix (a `p` by `p` matrix).
#' }
#' 
#' 
#' @examples
#' p = 500
#' n = 30
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' 
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' estimatedCov_sample = cov(X)
#' estimatedCov_analytical = cov_analytical_NL_shrinkage(X)
#' estimatedCov_QIS = cov_quadratic_inverse_shrinkage(X)
#' 
#' # We now compare the distances between the truth and the estimators.
#' LossEuclideanEigenvalues2(estimatedCov_sample, Sigma, type = "covariance")
#' LossEuclideanEigenvalues2(estimatedCov_analytical, Sigma)
#' LossEuclideanEigenvalues2(estimatedCov_QIS, Sigma)
#'
#' LossFrobenius2(estimatedCov_sample, Sigma, type = "covariance")
#' LossFrobenius2(estimatedCov_analytical, Sigma)
#' LossFrobenius2(estimatedCov_QIS, Sigma)
#' 
#' for (m in 1:4){
#'   estimatedCov_shrink_higher = cov_higher_order_shrinkage(X, m = m)
#'   
#'   loss = LossFrobenius2(estimatedCov_shrink_higher, Sigma)
#'   cat("m = ", m, ", loss = ", loss, "\n")
#'   
#'   # cat("alpha = ", estimatedCov_shrink_higher$alpha)
#'   cat("\n")
#' }
#' 
#' estimator = sum(diag(estimatedCov_sample)) / p * diag(nrow = p)
#' LossFrobenius2(estimator, Sigma, type = "covariance")
#' 
#' oracle_Estimator = cov_NL_oracle(X, Sigma = Sigma)
#' LossFrobenius2(oracle_Estimator, Sigma)
#' 
#' 
#' \donttest{
#' # Example of spiked covariance matrix model
#' 
#' p1 <- 0.2*p
#' p2 <- 0.4*p
#' X0 <- matrix(rnorm(p * 10*p), nrow = p, ncol = 10*p)
#' H <- eigen(X0 %*% t(X0))$vectors
#' 
#' eigenvalues = c(rep(1, p1) , rep(3, p2) , rep(5, p-p1-p2-1), sqrt(p))
#' D <- diag(eigenvalues)
#' sqD <- diag(sqrt(eigenvalues))
#' 
#' Sigma <- H %*% D %*% t(H)
#' }
#' 
#' @export
cov_higher_order_shrinkage <- function(X, centeredCov = TRUE, m, verbose = 0){
  
  if (verbose > 0){
    cat("Starting `cov_higher_order_shrinkage`...\n")
  }
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  if (verbose > 0){
    cat("*  m = ", m, "\n")
  }
  
  list_power_S = list()
  power_S = Ip
  
  for (k in 1:(2 * m)){
    if (k == 1){
      power_S = S
    } else {
      power_S = power_S %*% S
    }
    list_power_S[[k]] <- power_S
  }
  
  estimatedM = compute_M_covariance(m = m, c_n = cn, p = p,
                                    S = S, verbose = verbose,
                                    list_power_S = list_power_S)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = estimatedM$invM %*% estimatedM$hm
  if (verbose > 0){
    cat("Optimal alpha: \n")
    print(alpha)
  }
  
  estimated_covariance_matrix = alpha[1] * Ip
  power_S = Ip
  
  for (k in 1:m){
    estimated_covariance_matrix = estimated_covariance_matrix + 
      alpha[k + 1] * list_power_S[[k]]
  }
  
  result = list(
    estimated_covariance_matrix = estimated_covariance_matrix,
    invM = estimatedM$invM,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v
  )
  
  class(result) <- c("EstimatedCovarianceMatrix")
  
  return (result)
}


compute_M_covariance <- function(m, c_n, p, S, verbose, list_power_S){
  q1 = tr(S) / p
  
  if (m == 0){
    
    result = list(invM = 1, hm = q1)
    return (result)
  }
  
  # Now m is at least 1
  
  all_tr = compute_trace_power_matrix(list_power_S = list_power_S,
                                      m = 2 * m, verbose = verbose - 1)
  
  u = compute_HigherOrderCov_hat_u(m = m + 1, c_n = c_n, p = p,
                                   verbose = verbose - 1, all_tr = all_tr)
  
  f = compute_HigherOrderCov_hat_f(u = u, m = m, p = p, verbose = verbose - 1,
                                   all_tr = all_tr)
  
  hm = c(q1, f)
  
  old_invM = matrix(1, ncol = 1, nrow = 1)
  for (m in 1:m){
    invM = matrix(ncol = m + 1, nrow = m + 1)
    
    m_tilde = matrix(all_tr[m:(2 *m - 1)] / p, ncol = 1)
    # Transposed version:
    m_tilde_t = t(m_tilde)
    
    vector_correction_term = old_invM %*% m_tilde
    
    xi_n_m = all_tr[2 * m] - m_tilde_t %*% vector_correction_term
    xi_n_m = as.numeric(xi_n_m)
    
    corr_term = vector_correction_term %*% t(vector_correction_term) / xi_n_m
    
    invM[1:m, 1:m] = old_invM + corr_term
    invM[1:m, m+1] = - vector_correction_term / xi_n_m
    invM[m+1, 1:m] = - t(vector_correction_term) / xi_n_m
    invM[m+1, m+1] = 1 / xi_n_m
    
    old_invM = invM
    if (verbose > 0){
      cat("m =", m, "M^{-1} =\n")
      print(invM)
    }
  }
  
  return (list(invM = invM, hm = hm))
}


# This returns a vector of size m, where the j-th element corresponds to hat f_j.
# This computes f_1, ...., f_m. Not f_0.
compute_HigherOrderCov_hat_f <- function(u, m, p, verbose, all_tr){
  q = rep(NA, m + 1)
  f = rep(NA, m)
  
  q[1] = all_tr[1] / p
  
  all_Bell_polynomials = matrix(ncol = m + 1, nrow = m + 1)
  
  for (j in 1:(m+1)){
    for (k in 1:j){
      
      all_Bell_polynomials[j, k] = 
        kStatistics::e_eBellPol(j, k, c(u[1:(j - k + 1)], rep(0, k - 1)) )
      
      if (verbose > 0){
        cat("This is term j =", j, ", k =", k,
            ", Bell polynomial =", all_Bell_polynomials[j, k], "\n")
      }
    }
  }
  
  for (j in 2:(m+1)){
    q[j] = all_tr[j] / p
    
    for (k in 1:(j - 1)){
      
      Bell_polynomial = all_Bell_polynomials[j, k]
      corr_term = (-1)^(k+j) * (factorial(k) / factorial(j) ) * q[k] * Bell_polynomial
      
      q[j] = q[j] - corr_term
      
      if (verbose > 0){
        cat("This is term j =", j, ", k =", k,
            ", corr_term =", corr_term, "\n")
      }
    }
  }
  
  for (j in 1:m){
    f[j] = 0
    
    for (k in 1:j){
      
      Bell_polynomial = all_Bell_polynomials[j, k]
      add_term = (-1)^(k+j) * (factorial(k) / factorial(j) ) * q[k+1] * Bell_polynomial
      
      f[j] = f[j] + add_term
      
      if (verbose > 0){
        cat("This is term j =", j, ", k =", k,
            ", add_term =", add_term, "\n")
      }
    }
  }
  
  return (f)
}

# This returns a vector of size m, where the j-th element corresponds to hat u_j.
compute_HigherOrderCov_hat_u <- function(m, c_n, p, verbose, all_tr){
  u = rep(NA, m)
  u[1] = 1
  
  for (j in 2:m){
    u[j] <- (-1)^(j - 1) * factorial(j) * (c_n / p) * all_tr[j - 1]
  }
  return (u)
}


#' Compute all traces until m
#' 
#' @returns a vector of size `m` containing the trace of S^j for j = 1, ..., m.
#' 
#' @noRd
compute_trace_power_matrix <- function(list_power_S, m, verbose){
  all_tr = rep(NA, m)
  
  for (j in 1:m){
    all_tr[j] <- tr(list_power_S[[j]])
  }
  return (all_tr)
}
