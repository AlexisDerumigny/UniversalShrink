
estimator_d0_thetaknown <- function(Ip, Sn, t, Theta){
  
  result = t * tr( solve(Sn + t * Ip) %*% Theta )
  return (result)
}

estimator_d1_thetaknown <- function(Ip, Sn, t, Theta, p, cn){
  
  iS_ridge = solve(Sn + t * Ip)
  
  numerator = t * tr( iS_ridge %*% iS_ridge %*% Theta ) - t * estimator_d0_thetaknown(Ip, Sn, t, Theta)
  denominator = cn * ((1 / p) * tr( iS_ridge %*% iS_ridge %*% Theta ) 
                      - t^(-2) * (cn - 1) / cn)
  
  result = numerator / denominator
  
  return (result)
}

estimator_q1 <- function(Sn, Theta){
  result = tr(Sn %*% Theta)
  return (result)
}

estimator_q2 <- function(Sn, Theta, p, cn){
  result = tr(Sn %*% Sn %*% Theta) - cn * (1/p) * tr(Sn) * tr(Sn %*% Theta)
  return (result)
}

#' Estimator of d0(t0, Sigma  /p)
estimator_d0_1p_Sigma <- function(t0, hat_v_t0, cn){
  
  result = 1 / (cn * hat_v_t0) - t0 / cn
  return (result)
}

#' Estimator of d0(t0, Sigma^2 / p)
estimator_d0_1p_Sigma2 <- function(t0, hat_v_t0, cn, Sn, verbose = verbose){
  first_term = (1 / hat_v_t0) * (tr(Sn) / p)
  second_term = (1 / hat_v_t0) * ( (1 / (cn * hat_v_t0)) - t0 / cn)
  
  result = (1 / hat_v_t0) * (tr(Sn) / p - (1 / (cn * hat_v_t0)) + t0 / cn)
  
  if (verbose){
    cat("Estimator of d0(t0, Sigma^2 / p) : \n")
    cat("*  first_term = ", first_term, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  stopifnot(first_term - second_term == result)
  
  return (result)
}

#' Estimator of d0(t0, Sigma^2 * Pi_0 / p)
estimator_d0_1p_Sigma2_Pi0 <- function(t0, hat_v_t0, cn, Pi0, Ip, Sn, verbose){
  first_term = (1 / hat_v_t0) * (1 / p) * tr(Sn %*% Pi0)
  
  d0_t0_1p_Pi0 = estimator_d0_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Pi0)
  second_term = (1 / hat_v_t0^2) * (tr(Pi0) / p - d0_t0_1p_Pi0)
  
  result = first_term - second_term
  
  if (verbose){
    cat("Estimator of d0(t0, Sigma^2 * Pi_0 / p) : \n")
    cat("*  first_term = ", first_term, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  return (result)
}

#' Estimator of d1(t0, Sigma^2 / p)
estimator_d1_1p_Sigma2 <- function(t0, hat_v_t0, cn, Pi0, Ip, Sn){
  first_term = (1 / hat_v_t0^2)
  
  d1_t0_1p_Ip = estimator_d1_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Ip / p,
                                        p = p, cn = cn)
  
  second_term = tr(Sn) / p  +  d1_t0_1p_Ip  - 2 / (cn * hat_v_t0)  +  2 * t0 / cn
  
  result = first_term * second_term
  return (result)
}


best_alphabeta_ridge_shrinkage <- function(t0, cn, Pi0, Ip, Sn, verbose = verbose){
  
  hat_v_t0 = estimator_vhat_derivative(t = t0, m = 0, Sn = Sn, Ip = Ip, cn = cn)
  hat_vprime_t0 = estimator_vhat_derivative(t = t0, m = 1, Sn = Sn, Ip = Ip, cn = cn)
  
  d0_1p_Sigma = estimator_d0_1p_Sigma(t0 = t0, hat_v_t0 = hat_v_t0, cn = cn)
  
  d0_1p_Sigma2 = estimator_d0_1p_Sigma2(t0 = t0, hat_v_t0 = hat_v_t0, cn = cn,
                                        Sn = Sn, verbose = verbose - 1)
  
  d0_1p_Sigma2_Pi0 = estimator_d0_1p_Sigma2_Pi0(t0 = t0, hat_v_t0 = hat_v_t0,
                                                cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn,
                                                verbose = verbose - 1)
  
  d1_1p_Sigma2 = estimator_d1_1p_Sigma2(t0 = t0, hat_v_t0 = hat_v_t0, cn = cn,
                                        Pi0 = Pi0, Ip = Ip, Sn = Sn)
  
  q1 = estimator_q1(Sn = Sn, Theta = Pi0 / p)
  q2 = estimator_q2(Sn = Sn, Theta = Pi0 %*% Pi0 / p, p = p, cn = cn)
  
  if (verbose){
    cat("Estimators: \n")
    cat("*  hat_v_t0 = ", hat_v_t0, "\n")
    cat("*  hat_vprime_t0 = ", hat_vprime_t0, "\n")
    cat("*  d0_1p_Sigma = ", d0_1p_Sigma, "\n")
    cat("*  d0_1p_Sigma2 = ", d0_1p_Sigma2, "\n")
    cat("*  d0_1p_Sigma2_Pi0 = ", d0_1p_Sigma2_Pi0, "\n")
    cat("*  d1_1p_Sigma2 = ", d1_1p_Sigma2, "\n")
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n")
    cat("\n")
  }
  
  # Computation of alpha =======================================================
  numerator_alpha = d0_1p_Sigma * q2 - d0_1p_Sigma2_Pi0 * q1
  
  denominator_alpha_term1 = (t0^(-1) * d0_1p_Sigma2 + hat_vprime_t0 * d1_1p_Sigma2) * q2
  denominator_alpha_term2 = t0^(-1) * d0_1p_Sigma2_Pi0^2
  
  denominator_alpha = denominator_alpha_term1 - denominator_alpha_term2
  
  alpha = numerator_alpha / denominator_alpha
  
  # Computation of beta ========================================================
  numerator_beta_term1 = (t0^(-1) * d0_1p_Sigma2 + hat_vprime_t0 * d1_1p_Sigma2) * q1
  numerator_beta_term2 = t0^(-1) * d0_1p_Sigma * d0_1p_Sigma2_Pi0
    
  numerator_beta = numerator_beta_term1 - numerator_beta_term2
  
  # Note: beta has the same denominator as alpha, so it can be directly reused.
  beta = numerator_beta / denominator_alpha
  
  if (verbose){
    cat("Optimal values: \n")
    cat("*  numerator_alpha = ", numerator_alpha, "\n")
    cat("*  numerator_beta = ", numerator_beta, "\n")
    cat("*  denominator = ", denominator_alpha, "\n")
    cat("*  alpha = ", alpha, "\n")
    cat("*  beta = ", beta, "\n")
    cat("\n")
  }
  
  result = list(alpha = alpha, beta = beta)
  
  return (result)
}



estimator_vhat_derivative <- function(t, m, Sn, Ip, cn){
  
  term1 = (-1)^m * factorial(m) * cn
  
  iS_ridge = solve(Sn + t * Ip)
  
  # We put this matrix to the power m+1
  iS_ridge_power_m1 = Ip
  for (i in 1:(m+1)){
    iS_ridge_power_m1 = iS_ridge_power_m1 %*% iS_ridge
  }
    
  term2 = tr(iS_ridge_power_m1) / p - t^(- (m+1) ) * (cn - 1) / cn
  
  result = term1 * term2
  
  return (result)
}


#' @rdname ridge_target_identity_optimal
#' @export
ridge_target_general_semioptimal <- function (Y, centeredCov, t, Pi0, verbose = 2){
  
  if (verbose){
    cat("Starting `ridge_target_general_semioptimal`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  if (verbose){
    cat("*  n = ", n, "\n")
    cat("*  p = ", p, "\n")
    cat("*  t = ", t, "\n")
  }
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    if (verbose){
      cat("*  centered case\n")
    }
    
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
    
    c_n = p / (n-1)
  } else {
    if (verbose){
      cat("*  non-centered case\n")
    }
    
    S <- Y %*% t(Y)/n
    
    # Inverse companion covariance
    iY <- solve(t(Y) %*% Y / n)
    
    # Moore-Penrose inverse
    iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
    
    c_n = p / n
  }
  
  if (verbose){
    cat("*  c_n = ", c_n, "\n\n")
  }
  
  iS_ridge <- solve(S + t * Ip)
  
  best_alphabeta = best_alphabeta_ridge_shrinkage(t0 = t, cn = c_n,
                                                  Pi0 = Pi0, Ip = Ip, Sn = S,
                                                  verbose = verbose)
  
  alpha <- best_alphabeta$alpha
  beta <- best_alphabeta$beta
  
  iS_ShRt1 <- alpha * iS_ridge + beta * Ip
  
  
  return(list(
    estimated_precision_matrix = iS_ShRt1,
    alpha_optimal = alpha,
    beta_optimal = beta
  ) )
}

