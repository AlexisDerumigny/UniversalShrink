

MPR_target_identity_optimal <- function (X, centeredCov = TRUE, verbose = 2,
                                         eps = 1/(10^6), upp = pi/2 - eps, 
                                         initialValue = 1.5){
  if (verbose > 0){
    cat("Starting `MPR_target_identity_optimal`...\n")
  }
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  
  hL2MPr <- function(u){
    loss = loss_L2_MPR_optimal(t = tan(u), S = S, cn = cn, p = p, Ip = Ip)
    return(loss)
  }
  
  hL2MPR_max <- stats::optim(par = initialValue, hL2MPr,lower = eps, upper = upp,
                             method = "L-BFGS-B", control = list(fnscale = -1))
  u_MPR <- hL2MPR_max$par
  t <- tan(u_MPR)
  
  if (verbose > 0){
    cat("*  optimal t =", t ,"\n\n")
  }
  
  iS_ridge <- solve(S + t * Ip)
  
  best_alphabeta = best_alphabeta_MPR_shrinkage_identity(
    p = p, t = t, cn = cn, S = S, iS_ridge = iS_ridge, verbose = verbose)
  
  alpha <- best_alphabeta$alpha
  beta <- best_alphabeta$beta
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  MPR_target_identity = alpha * MPR_estimator + beta * Ip
  
  result = list(
    estimated_precision_matrix = MPR_target_identity,
    t_optimal = t,
    alpha_optimal = alpha,
    beta_optimal = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}



loss_L2_MPR_optimal <- function(t, S, cn, p, Ip){
  
  r = (cn - 1) / cn
  
  iS_Rt<-solve(S+t*Ip)
  
  trS1_t<-sum(diag(iS_Rt))/p
  trS2_t<-sum(diag(iS_Rt%*%iS_Rt))/p
  trS3_t<-sum(diag(iS_Rt%*%iS_Rt%*%iS_Rt))/p
  trS4_t<-sum(diag(iS_Rt%*%iS_Rt%*%iS_Rt%*%iS_Rt))/p
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - cn * q1^2
  
  hvt<-cn*(trS1_t-r/t)
  hvprt<--cn*(trS2_t-r/t/t)
  hvprprt<-2*cn*(trS3_t-r/t/t/t)
  hvprprprt<- -6*cn*(trS4_t-r/t/t/t/t)
  
  ihvt<-1/hvt
  ihvt_2<-ihvt^2
  hvprt_2<-hvprt^2
  
  d0Sig_t<-ihvt/cn-t/cn
  d0Sig2_t<-ihvt*(q1-d0Sig_t)
  d1Sig_t<- (ihvt_2+1/hvprt)/cn
  d1Sig2_t<-ihvt*(d0Sig2_t-d1Sig_t)
  d2Sig2_t<-ihvt*(d1Sig2_t-(ihvt^3+hvprprt/(hvprt^3)/2)/cn)
  d3Sig2_t<-ihvt*(d2Sig2_t-(ihvt^4+hvprprt^2/(hvprt^5)/2-hvprprprt/(hvprt^4)/6)/cn)
  
  hgs2Sig2_t<- -(hvprt_2*d2Sig2_t-hvprprt*d1Sig2_t/2) + 
    t*(hvprprprt*d1Sig2_t/6-hvprt*hvprprt*d2Sig2_t+d3Sig2_t*hvprt^3)
  
  num_ShMPRt<- -d1Sig_t*q2+d1Sig2_t*q1
  den_ShMPRt<-hgs2Sig2_t*q2-hvprt_2*d1Sig2_t^2
  
  L2_ShMPRt<-hvprt_2*num_ShMPRt^2/den_ShMPRt/q2
  
  return(L2_ShMPRt)  
}



MPR_target_identity_semioptimal <- function (X, centeredCov = TRUE, t, verbose = 2){
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  iS_ridge <- solve(S + t * Ip)
  
  best_alphabeta = best_alphabeta_MPR_shrinkage_identity(
    p = p, t = t, cn = cn, S = S, iS_ridge = iS_ridge, verbose = verbose)
  
  alpha <- best_alphabeta$alpha
  beta <- best_alphabeta$beta
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  MPR_target_identity = alpha * MPR_estimator + beta * Ip
  
  result = list(
    estimated_precision_matrix = MPR_target_identity,
    t = t,
    alpha_optimal = alpha,
    beta_optimal = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


best_alphabeta_MPR_shrinkage_identity <- function(p, t, cn, S, iS_ridge, verbose)
{
  trS1_t<-sum(diag(iS_ridge))/p
  trS2_t<-sum(diag(iS_ridge%*%iS_ridge))/p
  trS3_t<-sum(diag(iS_ridge%*%iS_ridge%*%iS_ridge))/p
  trS4_t<-sum(diag(iS_ridge%*%iS_ridge%*%iS_ridge%*%iS_ridge))/p
  
  r = (cn - 1) / cn
  
  hvt<-cn*(trS1_t-r/t)
  hvprt<--cn*(trS2_t-r/t/t)
  hvprprt<-2*cn*(trS3_t-r/t/t/t)
  hvprprprt<- -6*cn*(trS4_t-r/t/t/t/t)
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - cn * q1^2
  
  ihvt<-1/hvt
  ihvt_2<-ihvt^2
  hvprt_2<-hvprt^2
  d0Sig_t<-ihvt/cn-t/cn
  d0Sig2_t<-ihvt*(q1-d0Sig_t)
  d1Sig_t<- (ihvt_2+1/hvprt)/cn
  d1Sig2_t<-ihvt*(d0Sig2_t-d1Sig_t)
  d2Sig2_t<-ihvt*(d1Sig2_t-(ihvt^3+hvprprt/(hvprt^3)/2)/cn)
  d3Sig2_t<-ihvt*(d2Sig2_t-(ihvt^4+hvprprt^2/(hvprt^5)/2-hvprprprt/(hvprt^4)/6)/cn)
  
  first_term_s2_Sigma2 = -(hvprt_2*d2Sig2_t-hvprprt*d1Sig2_t/2)
  
  second_term_s2_Sigma2_1 = hvprprprt*d1Sig2_t/6
  second_term_s2_Sigma2_2 = hvprt*hvprprt*d2Sig2_t
  second_term_s2_Sigma2_3 = d3Sig2_t*hvprt^3
  
  second_term_s2_Sigma2 = t * 
    (second_term_s2_Sigma2_1 - second_term_s2_Sigma2_2 + second_term_s2_Sigma2_3)
  
  hgs2Sig2_t = first_term_s2_Sigma2 + second_term_s2_Sigma2
  
  if (verbose > 0){
    cat("Estimator of MPR s2(t, Sigma^2) : \n")
    cat("*  d1_1p_Sigma2 = ", d1Sig2_t, "\n")
    cat("*  d2_1p_Sigma2 = ", d2Sig2_t, "\n")
    cat("*  d3_1p_Sigma2 = ", d3Sig2_t, "\n")
    cat("*  first_term = ", first_term_s2_Sigma2, "\n")
    cat("*  second_term = ", second_term_s2_Sigma2, "\n")
    cat("*    second_term_1 = ", second_term_s2_Sigma2_1, "\n")
    cat("*    second_term_2 = ", second_term_s2_Sigma2_2, "\n")
    cat("*    second_term_3 = ", second_term_s2_Sigma2_3, "\n")
    cat("*  result = ", hgs2Sig2_t, "\n\n")
  }
  
  if (verbose > 1){
    cat("Estimators: \n")
    cat("*  hat_v_t0 = ",       hvt,       "\n")
    cat("*  hat_vprime_t0 = ",  hvprt,     "\n")
    cat("*  hat_vsecond_t0 = ", hvprprt,   "\n")
    cat("*  hat_vthird_t0 = ",  hvprprprt, "\n")
    
    cat("*  d0_1p_Sigma = ", d0Sig_t, "\n")
    cat("*  d0_1p_Sigma2 = ", d0Sig2_t, "\n")
    cat("*  d1_1p_Sigma = ", d1Sig_t, "\n")
    cat("*  d1_1p_Sigma2 = ", d1Sig2_t, "\n")
    
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n")
    
    cat("*  s2_Sigma2 = ", hgs2Sig2_t, "\n")
    
    cat("\n")
  }
  
  numerator_alpha_term1 = hvprt * d1Sig_t * q2
  numerator_alpha_term2 = hvprt * d1Sig2_t * q1
  
  numerator_beta_term1 = hgs2Sig2_t * q1
  numerator_beta_term2 = hvprt_2 * d1Sig_t * d1Sig2_t
  
  denominator_term1 = hgs2Sig2_t * q2
  denominator_term2 = hvprt_2 * d1Sig2_t^2
  
  numerator_alpha = - numerator_alpha_term1 + numerator_alpha_term2
  numerator_beta  =   numerator_beta_term1  - numerator_beta_term2
  denominator     = denominator_term1 - denominator_term2
  
  alpha = numerator_alpha / denominator
  beta  = numerator_beta  / denominator
  
  if (verbose > 0){
    cat("Optimal values: \n")
    
    cat("*  numerator_alpha = ", numerator_alpha, "\n")
    if (verbose > 1){
      cat("   *  first_term = ", numerator_alpha_term1, "\n")
      cat("   *  second_term = ", numerator_alpha_term2, "\n")
    }
    
    cat("*  numerator_beta = ", numerator_beta, "\n")
    if (verbose > 1){
      cat("   *  first_term = ", numerator_beta_term1, "\n")
      cat("   *  second_term = ", numerator_beta_term2, "\n")
    }
    
    cat("*  denominator = ", denominator, "\n")
    if (verbose > 1){
      cat("   *  first_term = ", denominator_term1, "\n")
      cat("   *  second_term = ", denominator_term2, "\n")
    }
    
    cat("*  alpha = ", alpha, "\n")
    cat("*  beta = ", beta, "\n")
    cat("\n")
  }
  
  result = list(alpha = alpha, beta = beta)
  
  return (result)
}


MPR_target_identity <- function (X, centeredCov = TRUE, t, alpha, beta, verbose = 0){
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  iS_ridge <- solve(S + t * Ip)
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  # This can also be written as:
  # iS_ridge %*% S %*% iS_ridge
  
  MPR_target_identity = alpha * MPR_estimator + beta * Ip
  
  result = list(
    estimated_precision_matrix = MPR_target_identity,
    t = t,
    alpha = alpha,
    beta = beta
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


