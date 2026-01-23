
#' Moore-Penrose-Ridge with identity target
#' 
#' This function computes
#' \deqn{\widehat{\Sigma^{-1}}^{ridge}_t
#'       = (S + t I_p)^{-1} - t * (S + t I_p)^{-2}},
#' where \eqn{S} is the sample covariance matrix and \eqn{t} is a given
#' parameter.
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
#' 
#' precision_MPR_Cent_optimal = MPR_target_identity_optimal(Y = Y,
#'                                                          centeredCov = TRUE)
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent_optimal, Sigma = Sigma),
#'     ", t opt = ", precision_MPR_Cent_optimal$t_optimal, 
#'     ", alpha opt = ", precision_MPR_Cent_optimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_Cent_optimal$beta_optimal, "\n", sep = "")
#' 
#' precision_MPR_Cent = MPR_target_identity(
#'    Y = Y, centeredCov = TRUE,
#'    t = precision_MPR_Cent_optimal$t_optimal,
#'    alpha = precision_MPR_Cent_optimal$alpha_optimal,
#'    beta = precision_MPR_Cent_optimal$beta_optimal)
#'    
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent, Sigma = Sigma))
#' 
#' precision_MPR_Cent = MPR_target_identity(
#'    Y = Y, centeredCov = TRUE,
#'    t = precision_MPR_Cent_optimal$t_optimal,
#'    alpha = 1, beta = 0)
#'    
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent, Sigma = Sigma))
#' 
#' precision_MPR_Cent = MPR_no_shrinkage(Y = Y, centeredCov = TRUE,
#'                                       t = precision_MPR_Cent_optimal$t_optimal)
#' cat("loss = ", FrobeniusLoss2(precision_MPR_Cent, Sigma = Sigma))
#' 
#' 
#' @export
MPR_target_identity_optimal <- function (Y, centeredCov = TRUE, verbose = 2,
                                         eps = 1/(10^6), upp = pi/2 - eps){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  
  hL2MPr <- function(u){
    loss = loss_L2_MPR_optimal(t = tan(u), S = S, cn = cn, p = p, Ip = Ip)
    return(loss)
  }
  
  
  hL2MPR_max <- stats::optim(1.5, hL2MPr,lower = eps, upper = upp,
                             method = "L-BFGS-B", control = list(fnscale = -1))
  u_MPR<- hL2MPR_max$par
  t <- tan(u_MPR)
  
  iS_ridge <- solve(S + t * Ip)
  
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
  hgs2Sig2_t<- -(hvprt_2*d2Sig2_t-hvprprt*d1Sig2_t/2) +
    t*(hvprprprt*d1Sig2_t/6-hvprt*hvprprt*d2Sig2_t+d3Sig2_t*hvprt^3)
  
  num_a_ShMPRt<- -hvprt*d1Sig_t*q2+hvprt*d1Sig2_t*q1
  num_b_ShMPRt<- hgs2Sig2_t*q1- hvprt_2*d1Sig_t*d1Sig2_t
  den_ShMPRt<-hgs2Sig2_t*q2-hvprt_2*d1Sig2_t^2
  alpha <- num_a_ShMPRt/den_ShMPRt
  beta  <- num_b_ShMPRt/den_ShMPRt
  
  
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


#' @rdname MPR_target_identity_optimal
#' @export
MPR_target_identity_semioptimal <- function (Y, centeredCov = TRUE, t, verbose = 2){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  cn = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  iS_ridge <- solve(S + t * Ip)
  
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
    cat("*  s2_Sigma2 = ", hgs2Sig2_t, "\n")
    cat("*    first_term = ", first_term_s2_Sigma2, "\n")
    cat("*    second_term = ", second_term_s2_Sigma2, "\n")
    cat("*      second_term_1 = ", second_term_s2_Sigma2_1, "\n")
    cat("*      second_term_2 = ", second_term_s2_Sigma2_2, "\n")
    cat("*      second_term_3 = ", second_term_s2_Sigma2_3, "\n")
  }
  
  numerator_alpha_1 = hvprt*d1Sig_t*q2
  numerator_alpha_2 = hvprt*d1Sig2_t*q1
  numerator_alpha   = - numerator_alpha_1 + numerator_alpha_2
  
  numerator_beta_1 = hgs2Sig2_t*q1
  numerator_beta_2 = hvprt_2*d1Sig_t*d1Sig2_t
  numerator_beta   = numerator_beta_1 - numerator_beta_2
  
  denominator_1 = hgs2Sig2_t*q2
  denominator_2 = hvprt_2*d1Sig2_t^2
  denominator <- denominator_1 - denominator_2
  
  alpha = numerator_alpha / denominator
  beta  = numerator_beta  / denominator
  
  if (verbose > 0){
    cat("Optimal values: \n")
    cat("*  numerator_alpha = ", numerator_alpha, "\n")
    cat("*    numerator_alpha_1 = ", numerator_alpha_1, "\n")
    cat("*    numerator_alpha_2 = ", numerator_alpha_2, "\n")
    cat("*  numerator_beta = ", numerator_beta, "\n")
    cat("*    numerator_beta_1 = ", numerator_beta_1, "\n")
    cat("*    numerator_beta_2 = ", numerator_beta_2, "\n")
    cat("*  denominator = ", denominator, "\n")
    cat("*    denominator_1 = ", denominator_1, "\n")
    cat("*    denominator_2 = ", denominator_2, "\n")
    cat("*  alpha = ", alpha, "\n")
    cat("*  beta = ", beta, "\n")
    cat("\n")
  }
  
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


#' @rdname MPR_target_identity_optimal
#' @export
MPR_target_identity <- function (Y, centeredCov = TRUE, t, alpha, beta, verbose = 0){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
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


