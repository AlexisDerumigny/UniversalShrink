



ridge_target_identity_optimal <- function (Y, centeredCov = TRUE, verbose = 0,
                                           eps = 1/(10^6), upp = pi/2 - eps,
                                           initialValue = 1.5){
  if (verbose > 0){
    cat("Starting `ridge_target_identity_optimal`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  
  r = (c_n - 1)/c_n
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  ##### shrinkage Ridge
  
  hL2R <- function(u)
  {
    t <- tan(u)
    iS_Rt <- solve(S + t * Ip)
    trS1_t <- sum(diag(iS_Rt)) / p
    trS2_t <- sum(diag(iS_Rt %*% iS_Rt)) / p
    
    hvt <- c_n * (trS1_t - r / t)
    hvprt <- - c_n * (trS2_t - r / t^2)
    ihvt <- 1 / hvt
    ihvt_2 <- ihvt^2
    
    d0_t <- t * trS1_t
    d1_t <- - ( t * trS2_t - d0_t / t) / (c_n * (trS2_t - r / t^2) )
    
    d0Sig_t <- ihvt / c_n - t / c_n
    d0Sig2_t <- ihvt * (q1 - d0Sig_t)
    d1Sig2_t <- ihvt_2 * (q1 + d1_t - 2 * d0Sig_t)
    
    num_a_ShRt <- d0Sig_t * q2 - d0Sig2_t * q1
    den_ShRt <- (d0Sig2_t + t * hvprt * d1Sig2_t) * q2 - d0Sig2_t^2
    L2_ShRt <- num_a_ShRt^2 / den_ShRt / q2
    
    return(L2_ShRt)
  }
  
  
  hL2R_max <- stats::optim(par = initialValue, fn = hL2R,
                           lower = eps, upper = upp,
                           method= "L-BFGS-B", control = list(fnscale = -1))
  
  u_R <- hL2R_max$par
  
  t_R <- tan(u_R)
  iS_Rt1 <- solve(S+t_R*Ip)
  trS1_t1 <- tr(iS_Rt1) / p
  trS2_t1 <- tr(iS_Rt1 %*% iS_Rt1) / p
  
  hvt1 <- c_n * (trS1_t1 - r / t_R)
  hvprt1 <- - c_n * (trS2_t1 - r / t_R^2)
  ihvt1 <- 1 / hvt1
  ihvt1_2 <- ihvt1^2
  
  d0_t1 <- t_R * trS1_t1
  d1_t1 <- -(t_R * trS2_t1 - d0_t1 / t_R) / ( c_n * (trS2_t1 - r / t_R^2) )
  
  d0Sig_t1 <- ihvt1 / c_n - t_R / c_n
  d0Sig2_t1 <- ihvt1 * (q1 - d0Sig_t1)
  d1Sig2_t1 <- ihvt1_2*(q1 + d1_t1 - 2 * d0Sig_t1)
  
  num_a_ShRt1 <- d0Sig_t1 * q2 - d0Sig2_t1 * q1
  num_b_ShRt1 <- (d0Sig2_t1 / t_R + hvprt1 * d1Sig2_t1) * q1 - d0Sig_t1 * d0Sig2_t1 / t_R
  den_ShRt1 <- (d0Sig2_t1 / t_R + hvprt1 * d1Sig2_t1) * q2 - d0Sig2_t1^2 / t_R
  ha_ShRt1 <- num_a_ShRt1 / den_ShRt1
  hb_ShRt1 <- num_b_ShRt1 / den_ShRt1
  iS_ShRt1 <- ha_ShRt1 * iS_Rt1 + hb_ShRt1 * Ip
  
  result = list(
    estimated_precision_matrix = iS_ShRt1,
    alpha_optimal = ha_ShRt1,
    beta_optimal = hb_ShRt1,
    t_optimal = t_R
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}



ridge_target_identity_semioptimal <- function (Y, centeredCov, t, verbose = 2){
  
  if (verbose > 0){
    cat("Starting `ridge_target_identity_semioptimal`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  if (verbose > 0){
    cat("*  t = ", t, "\n")
  }
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  
  r = (c_n - 1)/c_n
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  
  iS_ridge <- solve(S + t * Ip)
  trS1_t <- tr(iS_ridge) / p
  trS2_t <- tr(iS_ridge %*% iS_ridge) / p
  
  hvt1 <- c_n * (trS1_t - r / t)
  hvprt1 <- - c_n * (trS2_t - r / t^2)
  ihvt1 <- 1 / hvt1
  ihvt1_2 <- ihvt1^2
  
  d0_t1 <- t * trS1_t
  d1_t1 <- -(t * trS2_t - d0_t1 / t) / ( c_n * (trS2_t - r / t^2) )
  
  d0Sig_t1 <- ihvt1 / c_n - t / c_n
  d0Sig2_t1 <- ihvt1 * (q1 - d0Sig_t1)
  d1Sig2_t1 <- ihvt1_2*(q1 + d1_t1 - 2 * d0Sig_t1)
  
  num_a_ShRt1 <- d0Sig_t1 * q2 - d0Sig2_t1 * q1
  num_b_ShRt1 <- (d0Sig2_t1 / t + hvprt1 * d1Sig2_t1) * q1 - d0Sig_t1 * d0Sig2_t1 / t
  
  
  if (verbose > 0){
    cat("Estimators: \n")
    cat("*  hat_v_t0 = ", hvt1, "\n")
    cat("*  hat_vprime_t0 = ", hvprt1, "\n")
    cat("*  d0(t, Theta) = ", d0_t1, "\n")
    cat("*  d1(t, Theta) = ", d1_t1, "\n")
    cat("*  d0_1p_Sigma = ", d0Sig_t1, "\n")
    cat("*  d0_1p_Sigma2 = ", d0Sig2_t1, "\n")
    cat("*  d1_1p_Sigma2 = ", d1Sig2_t1, "\n")
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n")
    cat("\n")
  }
  
  den_ShRt1 <- (d0Sig2_t1 / t + hvprt1 * d1Sig2_t1) * q2 - d0Sig2_t1^2 / t
  
  # Equivalent representation in terms of hm and M
  hm = c(q1, d0Sig_t1 / t)
  
  M = c(q2            , d0Sig2_t1 / t,
        d0Sig2_t1 / t , d0Sig2_t1 / t^2 + hvprt1 * d1Sig2_t1 / t)
  
  M = matrix(M, ncol = 2, byrow = TRUE)
  
  if (verbose > 0){
    cat("Term of M[2,2] computed as:\n")
    cat("*  d0Sig2_t1 / t^2: ", d0Sig2_t1 / t^2, "\n")
    cat("*  hvprt1 * d1Sig2_t1 / t: ", hvprt1 * d1Sig2_t1 / t, "\n\n")
  }
  
  
  # det_M = den_ShRt1 / t
  # 
  # inv_M = matrix( 
  #   c(d0Sig2_t1 / t^2 + hvprt1 * d1Sig2_t1 / t, - d0Sig2_t1 / t,
  #     - d0Sig2_t1 / t                         , q2 ),
  #   ncol = 2, byrow = TRUE) / 
  #   det_M
  
  
  alpha <- num_a_ShRt1 / den_ShRt1
  beta <- num_b_ShRt1 / den_ShRt1
  
  iS_ShRt1 <- alpha * iS_ridge + beta * Ip
  
  if (verbose > 0){
    cat("Optimal values: \n")
    cat("*  numerator_alpha = ", num_a_ShRt1, "\n")
    cat("*  numerator_beta = ", num_b_ShRt1, "\n")
    cat("*  denominator = ", den_ShRt1, "\n")
    cat("*  alpha = ", alpha, "\n")
    cat("*  beta = ", beta, "\n")
    cat("\n")
  }
  
  result = list(
    estimated_precision_matrix = iS_ShRt1,
    alpha_optimal = alpha,
    beta_optimal = beta,
    M = M,
    hm = hm
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}



ridge_target_identity <- function (Y, centeredCov = TRUE, t, alpha, beta,
                                   verbose = 0){
  if (verbose > 0){
    cat("Starting `ridge_target_identity`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  if (verbose > 0){
    cat("*  t = ", t, "\n")
  }
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  iS_ridge <- solve(S + t * Ip)
  
  iS_ShRt1 <- alpha * iS_ridge + beta * Ip
  
  result = list(
    estimated_precision_matrix = iS_ShRt1
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

