
# Compute the matrix M for the higher-order shrinkage of the Moore-Penrose-Ridge
# estimator
# 
# @return a list with
# - a square matrix M of size (m + 1)
# - a vector hm of size (m + 1)
# - the estimator of v
compute_M_t_MPR <- function(m, c_n, S_t_inverse, q1, q2, t, method_invM, verbose,
                            mpfr, precBits)
{
  s_and_v = compute_sv_ridge(m = 2 * m, # we need additional powers compared
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
    for (k in 0:j){
      additional_term = (-t)^(j - k) * choose(n = j, k = k) * s_ridge[2 * j - k, ]
      s[j, ] = s[j, ] + additional_term
    }
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
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  if (method_invM == "solve"){
    alpha = solve(M) %*% hm
  } else if (method_invM == "ginv"){
    if (! requireNamespace("MASS", quietly = TRUE)){
      stop("MASS needs to be installed to use `method_invM == 'ginv'.`")
    }
    alpha = MASS::ginv(M) %*% hm
  } else if (method_invM == "recursive"){
    if (verbose > 0){
      cat("Using the recursive formula to compute the inverse of the matrix M...\n")
    }
    
    # We avoid computing M and inverting it numerically. Here we compute the
    # inverse of the matrix M by using the recursive formula.
    invM = compute_M_inverse(m = m, all_tr0 = 1 / s2[1],
                             all_tr = s2[-1], verbose = verbose - 2,
                             mpfr = mpfr, precBits = precBits)
    if (verbose > 1){
      cat("M^{-1} = \n")
      print(invM)
    }
    
    alpha = invM %*% hm
  } else {
    stop("method_invM '", method_invM, "' unavailable. Possible choices are: ",
         "'solve' and 'recursive'.")
  }
  
  if (verbose > 0){
    cat("Optimal alpha: \n")
    print(alpha)
  }
  
  
  return(list(M = M, hm = hm, v = v, alpha = alpha))
}



#' Moore-Penrose-Ridge higher order shrinkage
#' 
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' @param t \code{t} is the penalization parameter.
#' 
#' @inheritParams cov_with_centering
#' 
#' @template param-mpfr
#' @template param-optimizationControls
#' @template param-method_invM
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
#' precision_higher_order_shrinkage_Cent = MPR_higher_order_shrinkage(X, m = 2, verbose = 1)
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
MPR_higher_order_shrinkage <- function(
    X, m = 3, centeredCov = TRUE, t = NULL, optimizationControls = NULL,
    method_invM = "recursive", verbose = 0, mpfr = FALSE, precBits = 2^16)
{
  call_ = match.call()
  if (is.null(t)){
    result = MPR_higher_order_shrinkage_optimal(
      X = X, m = m, centeredCov = centeredCov, verbose = verbose,
      optimizationControls = optimizationControls, method_invM = method_invM,
      call_ = call_,
      mpfr = mpfr, precBits = precBits)
    
  } else {
    result = MPR_higher_order_shrinkage_non_optimized(
      X = X, m = m, centeredCov = centeredCov, t = t, verbose = verbose, 
      method_invM = method_invM, call_ = call_,
      mpfr = mpfr, precBits = precBits)
  }
}


MPR_higher_order_shrinkage_non_optimized <- function(
    X, m, centeredCov = TRUE, t, method_invM = "recursive", verbose = 0,
    call_ = NULL, mpfr = FALSE, precBits = 2^16)
{
  if (verbose > 0){
    cat("Starting `MPR_higher_order_shrinkage`...\n")
  }
  check_Rmpfr(mpfr)
  
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
  
  estimatedM = compute_M_t_MPR(
    m = m, c_n = c_n, q1 = q1, q2 = q2, S_t_inverse = iS_ridge,
    t = t, method_invM = method_invM, verbose = verbose,
    mpfr = mpfr, precBits = precBits)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  # TODO: loss = 1 - t(estimatedM$hm) %*% solve(estimatedM$M) %*% estimatedM$hm
  # as in ridge_higher_order_shrinkage_optimal
  
  alpha = estimatedM$alpha
  
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
    v = estimatedM$v,
    t = t,
    m = m,
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}



MPR_higher_order_shrinkage_optimal <- function(
    X, m, centeredCov = TRUE, t, verbose = 0, optimizationControls = NULL,
    method_invM = "recursive", call_ = NULL, mpfr = FALSE, precBits = 2^16)
{
  if (verbose > 0){
    cat("Starting `MPR_higher_order_shrinkage_optimal` (with unknown t)...\n")
  }
  check_Rmpfr(mpfr)
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  if (verbose > 0){
    cat("*  m = ", m, "\n")
  }
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  
  estimatedLoss <- function(t){
    loss = loss_L2_MPR_higher_order_optimal(
      X = X, m = m, t = t, c_n = c_n, q1 = q1, q2 = q2, S = S, Ip = Ip,
      method_invM = method_invM, verbose = verbose - 1,
      mpfr = mpfr, precBits = precBits)
    
    return (loss)
  }
  
  if (is.null(optimizationControls)) {
    optimizationControls = list(method = "smoothed")
  }
  if (optimizationControls$method == "smoothed") {
    if (is.null(optimizationControls$grid)) {
      optimizationControls$grid <- grid_optimization_default(
        S = S, c_n = c_n, p = p, n = n)
    }
    
    if (is.null(optimizationControls$k)) {
      optimizationControls$k <- 2 * m + 1
    }
  }
  
  result_optimization = optimization(
    FUN = estimatedLoss, optimizationControls = optimizationControls,
    maximum = FALSE, verbose = verbose - 2)
  
  optimal_t = result_optimization$optimal_t
  
  if (verbose > 0){
    cat("*  optimal_t = ", optimal_t, "\n")
  }
  
  # ============================================================================
  # We now compute the estimator using this optimal_t that was found.
  
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  iS_ridge <- solve(S + optimal_t * Ip)
  
  MPR_estimator <- iS_ridge - optimal_t * iS_ridge %*% iS_ridge
  
  estimatedM = compute_M_t_MPR(
    m = m, c_n = c_n, q1 = q1, q2 = q2, S_t_inverse = iS_ridge,
    t = optimal_t, method_invM = method_invM, verbose = verbose - 2,
    mpfr = mpfr, precBits = precBits)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = estimatedM$alpha
  
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
    v = estimatedM$v,
    m = m,
    optimal_t = optimal_t,
    result_optimization = result_optimization,
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


loss_L2_MPR_higher_order_optimal <- function (
    X, m, t, c_n, q1, q2, S, Ip, method_invM = "recursive", verbose, 
    mpfr, precBits)
{
  if (verbose > 0){
    cat("t = ", t)
    if (verbose > 1){
      cat("\n")
    } else {
      cat(", ")
    }
  }
  
  S_t_inverse = solve(S + t * Ip)
    
  loss = tryCatch({
    estimatedM = compute_M_t_MPR(
      m = m, c_n = c_n, q1 = q1, q2 = q2, S_t_inverse = S_t_inverse,
      t = t, method_invM = method_invM, verbose = verbose - 2,
      mpfr = mpfr, precBits = precBits)
    
    loss = 1 - t(estimatedM$hm) %*% estimatedM$alpha
  }, error = function(e){e}
  )
  
  if (inherits(loss, "simpleError") || !is.finite(loss)){
    loss <- .Machine$double.xmax
  }
  
  if (verbose > 0){
    cat("loss = ", loss, "\n")
  }
  
  return (loss)
}

