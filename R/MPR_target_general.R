

#' Moore-Penrose-Ridge with general target
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
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the precision matrix, of class
#' `EstimatedPrecisionMatrix`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2024).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
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
#' 
#' t = 1
#' 
#' precision_MPR_Cent_gen_semioptimal = 
#'   MPR_target_general_semioptimal(Y = Y, centeredCov = TRUE, Pi0 = Ip, t = t)
#'                                                                 
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent_gen_semioptimal, Sigma = Sigma),
#'     ", t = ", t, 
#'     ", alpha opt = ", precision_MPR_Cent_gen_semioptimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_Cent_gen_semioptimal$beta_optimal, "\n", sep = "")
#' 
#' precision_MPR_Cent_semioptimal = 
#'   MPR_target_identity_semioptimal(Y = Y, centeredCov = TRUE, t = t)
#'   
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent_semioptimal, Sigma = Sigma),
#'     ", t = ", t, 
#'     ", alpha opt = ", precision_MPR_Cent_semioptimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_Cent_semioptimal$beta_optimal, "\n", sep = "")
#' 
#' precision_MPR_Cent_gen_optimal = 
#'   MPR_target_general_optimal(Y = Y, centeredCov = TRUE, Pi0 = Ip)
#'   
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent_gen_optimal, Sigma = Sigma),
#'     ", t = ", precision_MPR_Cent_gen_optimal$t_optimal, 
#'     ", alpha opt = ", precision_MPR_Cent_gen_optimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_Cent_gen_optimal$beta_optimal, "\n", sep = "")
#' 
#' precision_MPR_Cent_gen_optimal_oracle = 
#'   MPR_target_general_optimal(Y = Y, centeredCov = TRUE,
#'   Pi0 = solve(0.99 * Sigma + 0.01 * Ip))
#'   
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent_gen_optimal_oracle, Sigma = Sigma),
#'     ", t = ", precision_MPR_Cent_gen_optimal_oracle$t_optimal, 
#'     ", alpha opt = ", precision_MPR_Cent_gen_optimal_oracle$alpha_optimal,
#'     ", beta opt = ", precision_MPR_Cent_gen_optimal_oracle$beta_optimal, "\n", sep = "")
#' 
#' @export
#' 
MPR_target_general_optimal <- function (Y, centeredCov, Pi0, verbose = 3){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    cn = p / (n-1)
  } else {
    S <- Y %*% t(Y)/n
    
    cn = p / n
  }
  
  
  # TODO: provide this as an option for the user
  initialValue = 1.5
  eps <- 1/(10^6)
  upp <- pi/2 - eps
  
  hL2R <- function(u){
    loss = loss_L2_MPR_optimal_target_general(
      t = tan(u), Sn = S, p = p, Ip = Ip, cn = cn, Pi0 = Pi0,
      verbose = verbose - 3)
    
    if (verbose > 1){
      cat("t = ", tan(u), "; loss = ", loss, ". ")
    }
    
    return(loss)
  }
  
  control = list(fnscale = -1,
                 trace = if(verbose > 2){6} else {0},
                 factr = 1e8
                 # ndeps = 0.01
                 )
  
  hL2R_max <- stats::optim(par = initialValue, fn = hL2R,
                           lower = eps, upper = upp,
                           method = "L-BFGS-B", control = control)
  
  u_R <- hL2R_max$par
  t <- tan(u_R)
  if (verbose > 0){
    cat("*  optimal t =", t ,"\n")
  }
  
  iS_ridge <- solve(S + t * Ip)
  
  best_alphabeta = best_alphabeta_MPR_shrinkage(p = p, t0 = t, cn = cn, Pi0 = Pi0,
                                                Ip = Ip, Sn = S, verbose = verbose)
  
  alpha <- best_alphabeta$alpha
  beta <- best_alphabeta$beta
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  MPR_target_general = alpha * MPR_estimator + beta * Pi0
  
  result = list(
    estimated_precision_matrix = MPR_target_general,
    t_optimal = t,
    alpha_optimal = alpha,
    beta_optimal = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}




#' @param t the value of that parameter
#' @param Sn the sample covariance matrix (potentially centered)
#' @param Ip the identity matrix
#' @param cn the ratio p/n (potentially centered)
#' @param Pi0 the target
#'
#' @returns an estimator of the L2 loss
#'
#' @noRd
loss_L2_MPR_optimal_target_general <- function(t, Sn, p, Ip, cn, Pi0, verbose)
{
  hat_v_t0 = estimator_vhat_derivative(t = t, m = 0, Sn = Sn, p = p, 
                                       Ip = Ip, cn = cn)
  hat_vprime_t0 = estimator_vhat_derivative(t = t, m = 1, Sn = Sn, p = p, 
                                            Ip = Ip, cn = cn)
  hat_vsecond_t0 = estimator_vhat_derivative(t = t, m = 2, Sn = Sn, p = p, 
                                             Ip = Ip, cn = cn)
  hat_vthird_t0 = estimator_vhat_derivative(t = t, m = 3, Sn = Sn, p = p, 
                                            Ip = Ip, cn = cn)
  
  q1 = estimator_q1(Sn = Sn, Theta = Pi0 / p)
  q2 = estimator_q2(Sn = Sn, Theta = Pi0 %*% Pi0 / p, p = p, cn = cn)
  
  d1_1p_Sigma = estimator_d1_1p_Sigma(hat_v_t0 = hat_v_t0,
                                      hat_vprime_t0 = hat_vprime_t0, cn = cn)
  
  d1_1p_Sigma2_Pi0 = estimator_d1_1p_Sigma2Pi0(t0 = t, hat_v_t0 = hat_v_t0,
                                               cn = cn, p = p, Ip = Ip, Sn = Sn,
                                               Pi0 = Pi0, verbose = verbose)
  
  s2_Sigma2 = estimator_MPR_s2_Sigma2(hat_v_t0 = hat_v_t0,
                                      hat_vprime_t0 = hat_vprime_t0,
                                      hat_vsecond_t0 = hat_vsecond_t0,
                                      hat_vthird_t0 = hat_vthird_t0,
                                      t0 = t, p = p, cn = cn, Pi0 = Pi0, Ip = Ip,
                                      Sn = Sn, verbose = verbose)
  
  numerator = hat_vprime_t0^2 * (d1_1p_Sigma * q2 - d1_1p_Sigma2_Pi0 * q1)^2
  
  denominator_main = s2_Sigma2 * q2 - hat_vprime_t0^2 * d1_1p_Sigma2_Pi0^2
  
  denominator = q2 * denominator_main
  
  loss_L2 = numerator / denominator
  
  return (loss_L2)
}


#' @rdname MPR_target_general_optimal
#' @export
MPR_target_general_semioptimal <- function (Y, centeredCov, t, Pi0, verbose = 2){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    cn = p / (n-1)
  } else {
    S <- Y %*% t(Y)/n
    
    cn = p / n
  }
  
  
  iS_ridge <- solve(S + t * Ip)
  
  best_alphabeta = best_alphabeta_MPR_shrinkage(p = p, t0 = t, cn = cn, Pi0 = Pi0,
                                                Ip = Ip, Sn = S, verbose = verbose)
  
  alpha <- best_alphabeta$alpha
  beta <- best_alphabeta$beta
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  MPR_target_general = alpha * MPR_estimator + beta * Pi0
  
  result = list(
    estimated_precision_matrix = MPR_target_general,
    t = t,
    alpha_optimal = alpha,
    beta_optimal = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


#' Estimator of d1(t0, Sigma / p)
#' @noRd
estimator_d1_1p_Sigma <- function(hat_v_t0, hat_vprime_t0, cn){
  first_term = 1 / hat_v_t0^2
  
  second_term = 1 / hat_vprime_t0
  
  result = (first_term + second_term) / cn
  return (result)
}

# Estimator of d1(t, Sigma^2 Pi0)
estimator_d1_1p_Sigma2Pi0 <- function(t0, hat_v_t0, cn, p, Ip, Sn, Pi0, verbose){
  
  d0_1p_Sigma2_Pi0 = estimator_d0_1p_Sigma2_Pi0(p = p, t0 = t0, hat_v_t0 = hat_v_t0,
                                                cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn, 
                                                verbose = verbose)
  
  d1_1p_Pi0 = estimator_ridge_d1_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Pi0 / p,
                                      p = p, cn = cn)
  
  d0_1p_Pi0 = estimator_ridge_d0_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Pi0 / p)
  
  first_term = d0_1p_Sigma2_Pi0 / hat_v_t0
  second_term = d1_1p_Pi0 / hat_v_t0^2
  third_term = ( (tr(Pi0) / p) - d0_1p_Pi0 ) / hat_v_t0^3
  
  result = first_term + second_term - third_term
  
  if (verbose > 0){
    cat("Estimator of d1(t, Sigma^2 Pi0) : \n")
    cat("*  first_term = ", first_term, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*  third_term = ", third_term, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  return (result)
}

# Estimator of d2(t, Sigma^2 / p)
estimator_d2_1p_Sigma2 <- function(hat_v_t0, hat_vprime_t0, hat_vsecond_t0,
                                   t0, p, cn, Pi0, Ip, Sn, verbose){
  
  d1_1p_Sigma2 = estimator_d1_1p_Sigma2(t0 = t0, hat_v_t0 = hat_v_t0, p = p,
                                        cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn)
  
  second_term_1 = 1 / hat_v_t0^3
  second_term_2 = hat_vsecond_t0 / (2 * hat_vprime_t0^3)
  second_term = (second_term_1 + second_term_2) / cn
  
  result = (d1_1p_Sigma2 - second_term) / hat_v_t0
  
  if (verbose > 0){
    cat("Estimator of d2(t, Sigma^2 / p) : \n")
    cat("*  first_term = ", d1_1p_Sigma2, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  return (result)
}

# Estimator of d3(t, Sigma^2 / p)
estimator_d3_1p_Sigma2 <- function(hat_v_t0, hat_vprime_t0, hat_vsecond_t0, hat_vthird_t0,
                                   t0, p, cn, Pi0, Ip, Sn, verbose){
  
  d2_1p_Sigma2 = estimator_d2_1p_Sigma2(hat_v_t0 = hat_v_t0,
                                        hat_vprime_t0 = hat_vprime_t0,
                                        hat_vsecond_t0 = hat_vsecond_t0,
                                        t0 = t0, p = p, cn = cn, Pi0 = Pi0,
                                        Ip = Ip, Sn = Sn, verbose = verbose - 1)
  
  second_term_1 = 1 / hat_v_t0^4
  second_term_2 = hat_vsecond_t0^2 / (2 * hat_vprime_t0^5)
  second_term_3 = hat_vthird_t0 / (6 * hat_vprime_t0^4)
  
  second_term = (second_term_1 + second_term_2 - second_term_3) / cn
  
  result = (d2_1p_Sigma2 - second_term) / hat_v_t0
  
  if (verbose > 0){
    cat("Estimator of d2(t, Sigma^2 / p) : \n")
    cat("*  first_term = ", d2_1p_Sigma2, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  return (result)
}


# Estimator of s2(t, Sigma^2)
estimator_MPR_s2_Sigma2 <- function(hat_v_t0, hat_vprime_t0, hat_vsecond_t0, hat_vthird_t0,
                                    t0, p, cn, Pi0, Ip, Sn, verbose){
  
  d1_1p_Sigma2 = estimator_d1_1p_Sigma2(t0 = t0, hat_v_t0 = hat_v_t0, p = p,
                                        cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn)
  
  d2_1p_Sigma2 = estimator_d2_1p_Sigma2(hat_v_t0 = hat_v_t0,
                                        hat_vprime_t0 = hat_vprime_t0,
                                        hat_vsecond_t0 = hat_vsecond_t0,
                                        t0 = t0, p = p, cn = cn, Pi0 = Pi0,
                                        Ip = Ip, Sn = Sn, verbose = verbose - 1)
  
  d3_1p_Sigma2 = estimator_d3_1p_Sigma2(hat_v_t0 = hat_v_t0,
                                        hat_vprime_t0 = hat_vprime_t0,
                                        hat_vsecond_t0 = hat_vsecond_t0,
                                        hat_vthird_t0 = hat_vthird_t0,
                                        t0 = t0, p = p, cn = cn, Pi0 = Pi0,
                                        Ip = Ip, Sn = Sn, verbose = verbose - 1)
  
  first_term_1 = hat_vprime_t0^2 * d2_1p_Sigma2
  first_term_2 = hat_vsecond_t0 * d1_1p_Sigma2 / 2
  first_term = - (first_term_1 - first_term_2)
  
  second_term_1 = hat_vthird_t0 * d1_1p_Sigma2 / 6
  second_term_2 = hat_vprime_t0 * hat_vsecond_t0 * d2_1p_Sigma2
  second_term_3 = hat_vprime_t0^3 * d3_1p_Sigma2
  second_term = t0 * (second_term_1 - second_term_2 + second_term_3)
  
  result = first_term + second_term
  
  if (verbose > 0){
    cat("Estimator of MPR s2(t, Sigma^2) : \n")
    cat("*  d1_1p_Sigma2 = ", d1_1p_Sigma2, "\n")
    cat("*  d2_1p_Sigma2 = ", d2_1p_Sigma2, "\n")
    cat("*  d3_1p_Sigma2 = ", d3_1p_Sigma2, "\n")
    cat("*  first_term = ", first_term, "\n")
    cat("*  second_term = ", second_term, "\n")
    cat("*    second_term_1 = ", second_term_1, "\n")
    cat("*    second_term_2 = ", second_term_2, "\n")
    cat("*    second_term_3 = ", second_term_3, "\n")
    cat("*  result = ", result, "\n\n")
  }
  
  return (result)
}


best_alphabeta_MPR_shrinkage <- function(p, t0, cn, Pi0, Ip, Sn, verbose = verbose){
  
  hat_v_t0 = estimator_vhat_derivative(t = t0, m = 0, Sn = Sn, p = p,
                                       Ip = Ip, cn = cn)
  hat_vprime_t0 = estimator_vhat_derivative(t = t0, m = 1, Sn = Sn, p = p,
                                            Ip = Ip, cn = cn)
  hat_vsecond_t0 = estimator_vhat_derivative(t = t0, m = 2, Sn = Sn, p = p,
                                             Ip = Ip, cn = cn)
  hat_vthird_t0 = estimator_vhat_derivative(t = t0, m = 3, Sn = Sn, p = p, 
                                            Ip = Ip, cn = cn)
  
  d0_1p_Sigma = estimator_d0_1p_Sigma(t0 = t0, hat_v_t0 = hat_v_t0, cn = cn)
  
  d0_1p_Sigma2 = estimator_d0_1p_Sigma2(p = p, t0 = t0, hat_v_t0 = hat_v_t0, cn = cn,
                                        Sn = Sn, verbose = verbose - 1)
  
  d0_1p_Sigma2_Pi0 = estimator_d0_1p_Sigma2_Pi0(p = p, t0 = t0, hat_v_t0 = hat_v_t0,
                                                cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn,
                                                verbose = verbose - 1)
  
  d1_1p_Sigma = estimator_d1_1p_Sigma(hat_v_t0 = hat_v_t0, hat_vprime_t0 = hat_vprime_t0,
                                      cn = cn)
  
  d1_1p_Sigma2 = estimator_d1_1p_Sigma2(t0 = t0, hat_v_t0 = hat_v_t0, p = p,
                                        cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn)
  
  d1_1p_Sigma2Pi0 = estimator_d1_1p_Sigma2Pi0(t0 = t0, hat_v_t0 = hat_v_t0,
                                              cn = cn, p = p, Ip = Ip, Sn = Sn,
                                              Pi0 = Pi0, verbose = verbose - 1)
  
  q1 = estimator_q1(Sn = Sn, Theta = Pi0 / p)
  q2 = estimator_q2(Sn = Sn, Theta = Pi0 %*% Pi0 / p, p = p, cn = cn)
  
  s2_Sigma2 = estimator_MPR_s2_Sigma2(hat_v_t0 = hat_v_t0, hat_vprime_t0 = hat_vprime_t0,
                                      hat_vsecond_t0 = hat_vsecond_t0,
                                      hat_vthird_t0 = hat_vthird_t0,
                                      t0 = t0, p = p, cn = cn, Pi0 = Pi0, Ip = Ip, Sn = Sn,
                                      verbose = verbose - 1)
  
  if (verbose > 0){
    cat("Estimators: \n")
    cat("*  hat_v_t0 = ", hat_v_t0, "\n")
    cat("*  hat_vprime_t0 = ", hat_vprime_t0, "\n")
    
    d0_t0_1p_Pi0 = estimator_ridge_d0_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Pi0 / p)
    
    d1_t0_1p_Ip = estimator_ridge_d1_thetaknown(Ip = Ip, Sn = Sn, t = t0, Theta = Ip / p,
                                          p = p, cn = cn)
    
    cat("*  d0(t, Theta) = ", d0_t0_1p_Pi0, "\n")
    cat("*  d1(t, Theta) = ", d1_t0_1p_Ip, "\n")
    
    cat("*  d0_1p_Sigma = ", d0_1p_Sigma, "\n")
    cat("*  d0_1p_Sigma2 = ", d0_1p_Sigma2, "\n")
    cat("*  d0_1p_Sigma2_Pi0 = ", d0_1p_Sigma2_Pi0, "\n")
    cat("*  d1_1p_Sigma = ", d1_1p_Sigma, "\n")
    cat("*  d1_1p_Sigma2 = ", d1_1p_Sigma2, "\n")
    cat("*  d1_1p_Sigma2_Pi0 = ", d1_1p_Sigma2Pi0, "\n")
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n")
    cat("*  s2_Sigma2 = ", s2_Sigma2, "\n")
    cat("\n")
  }
  
  # Computation of alpha =======================================================
  
  numerator_alpha_term1 = hat_vprime_t0 * d1_1p_Sigma * q2
  numerator_alpha_term2 = hat_vprime_t0 * d1_1p_Sigma2Pi0 * q1
  numerator_alpha = - numerator_alpha_term1 + numerator_alpha_term2
    
  denominator_alpha_term1 = s2_Sigma2 * q2
  denominator_alpha_term2 = hat_vprime_t0^2 * d1_1p_Sigma2Pi0^2
  
  denominator_alpha = denominator_alpha_term1 - denominator_alpha_term2
  
  alpha = numerator_alpha / denominator_alpha
  
  # Computation of beta ========================================================
  numerator_beta_term1 = s2_Sigma2 * q1
  numerator_beta_term2 = hat_vprime_t0^2 * d1_1p_Sigma * d1_1p_Sigma2Pi0
  
  numerator_beta = numerator_beta_term1 - numerator_beta_term2
  
  # Note: beta has the same denominator as alpha, so it can be directly reused.
  beta = numerator_beta / denominator_alpha
  
  if (verbose > 0){
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



#' @rdname MPR_target_general_optimal
#' @export
MPR_target_general <- function (Y, centeredCov, t, alpha, beta, Pi0){
  
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
  
  MPR_target_general = alpha * MPR_estimator + beta * Pi0
  
  result = list(
    estimated_precision_matrix = MPR_target_general,
    t = t,
    alpha = alpha,
    beta = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

