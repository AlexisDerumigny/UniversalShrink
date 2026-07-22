

#' Oracle estimator of the precision matrix
#' 
#' 
#' @param nameEstimator name of the estimator of which we compute the oracle of.
#' If \code{nameEstimator = "Moore-Penrose"}, this computes the finite-sample
#' oracle of the higher-order shrinked Moore-Penrose estimator. Other
#' possibilities are \code{"ridge"} and \code{"MPR"}.
#' 
#' @param Sigma true covariance matrix
#' 
#' @template param-mpfr
#' 
#' @inheritParams MPR_higher_order_shrinkage
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
#' for (m in 1:3){
#'   cat("\nm = ", m, "\n")
#'   cat("MP: \n")
#'   precision_higher_order_shrinkage_Cent = Moore_Penrose_higher_order_shrinkage(
#'       X, m = m, centeredCov = TRUE, verbose = 0)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "Moore-Penrose",
#'       centeredCov = TRUE,method_invM = "recursive", verbose = 0)
#'       
#'   print(LossInverseFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   print(LossInverseFrobenius2(oracle, Sigma))
#'   
#'   cat("ridge: \n")
#'   precision_higher_order_shrinkage_Cent = ridge_higher_order_shrinkage(
#'       X, m = m, centeredCov = TRUE, verbose = 0)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "ridge",
#'       centeredCov = TRUE,method_invM = "recursive", verbose = 0)
#'       
#'   print(LossInverseFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   print(LossInverseFrobenius2(oracle, Sigma))
#'   cat("optimal t = ", precision_higher_order_shrinkage_Cent$t,
#'       " (BF) ,  ", oracle$optimal_t, " (oracle) \n")
#'   
#'   cat("MPR: \n")
#'   precision_higher_order_shrinkage_Cent = MPR_higher_order_shrinkage(
#'       X, m = m, centeredCov = TRUE, verbose = 0)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "MPR",
#'       centeredCov = TRUE,method_invM = "recursive", verbose = 0)
#'       
#'   print(LossInverseFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   print(LossInverseFrobenius2(oracle, Sigma))
#'   cat("optimal t = ", precision_higher_order_shrinkage_Cent$optimal_t,
#'       " (BF) ,  ", oracle$optimal_t, " (oracle) \n")
#' }
#' 
#' 
#' @export
oracle_higher_order_shrinkage <- function(
    X, m, Sigma, nameEstimator = "Moore-Penrose", centeredCov = TRUE,
    method_invM = "recursive", verbose = 0, optimizationControls = NULL,
    mpfr = FALSE, precBits = 2^16) {
  
  call_ = match.call()
  check_Rmpfr(mpfr)
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (nameEstimator == "Moore-Penrose") {
    
    # Moore-Penrose inverse of the sample covariance matrix
    iS_MP <- Moore_Penrose(X = X, centeredCov = centeredCov)
    estimated_precision_matrix <- as.matrix(iS_MP)
    
  } else if (nameEstimator %in% c("ridge", "MPR")) {
    
    if (nameEstimator == "ridge") {
      
      estimatedLoss <- function(t){
        loss = loss_L2_ridge_oracle_higher_order_optimal(
          S = S, t = t, Ip = Ip, p = p, m = m, Sigma = Sigma,
          method_invM = method_invM, verbose = verbose - 2,
          mpfr = mpfr, precBits = precBits)
        
        return (loss)
      }
      
    } else if (nameEstimator == "MPR") {
      
      estimatedLoss <- function(t){
        loss = loss_L2_MPR_oracle_higher_order_optimal(
          S = S, t = t, Ip = Ip, p = p, m = m, Sigma = Sigma,
          method_invM = method_invM, verbose = verbose - 2,
          mpfr = mpfr, precBits = precBits)
        
        return (loss)
      }
    }
    
    if (is.null(optimizationControls)) {
      optimizationControls = list(method = "optimize")
    }
    if (optimizationControls$method == "smoothed" && 
        is.null(optimizationControls$grid) ) {
      
      optimizationControls$grid <- grid_optimization_default(S = S, c_n = c_n)
    }
    
    result_optimization = optimization(
      FUN = estimatedLoss, optimizationControls = optimizationControls,
      maximum = FALSE, verbose = verbose - 2)
    
    optimal_t = result_optimization$optimal_t
    
    if (verbose > 0){
      cat("*  optimal_t = ", optimal_t, "\n")
    }
    
    iS_ridge <- solve(S + optimal_t * Ip)
    
    if (nameEstimator == "ridge") {
      estimated_precision_matrix = iS_ridge
      
    } else if (nameEstimator == "MPR") {
      
      MPR_estimator <- iS_ridge - optimal_t * iS_ridge %*% iS_ridge
      estimated_precision_matrix = MPR_estimator
    }
    
  } else {
    stop("nameEstimator = '", nameEstimator, "' is not valid.",
         "It should be one of the following: 'Moore-Penrose', 'ridge', 'MPR'.")
  }
  
  
  resultM <- compute_M_oracle(
    estimator_S = estimated_precision_matrix, Sigma = Sigma, m = m, Ip = Ip,
    p = p, verbose = verbose, mpfr = mpfr, precBits = precBits)
  
  alpha = resultM$alpha
  
  result = alpha[1] * Ip
  power_estimated_precision = Ip
  
  for (k in 1:m){
    power_estimated_precision = 
      power_estimated_precision %*% estimated_precision_matrix
    
    result = result + alpha[k + 1] * power_estimated_precision
  }
  
  result = list(
    estimated_precision_matrix = result,
    M = resultM$M,
    hm = resultM$hm,
    # invM_solve = resultM$invM_solve,
    invM_recursive = resultM$invM_recursive,
    alpha = alpha,
    optimal_t = if(nameEstimator %in% c("ridge", "MPR")) {optimal_t} ,
    result_optimization = if(nameEstimator %in% c("ridge", "MPR")) {
      result_optimization} ,
    method = "Oracle higher-order shrinkage",
    nameBaselinEstimator = nameEstimator,
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

compute_M_oracle <- function(estimator_S, Sigma, m, Ip, p, verbose,
                             mpfr, precBits) {
  hm = rep(NA, m + 1)
  s2 = rep(NA, 2 * m + 1)
  
  power_estimator_S_Sigma = Sigma / p
  power_estimator_S_Sigma2 = Sigma %*% Sigma / p
  for (j in 1:(m + 1)){
    hm[j] = tr(power_estimator_S_Sigma)
    
    power_estimator_S_Sigma = estimator_S %*% power_estimator_S_Sigma
  }
  
  for (j in 1:(2 * m + 1)){
    s2[j] = tr(power_estimator_S_Sigma2)
    
    power_estimator_S_Sigma2 = estimator_S %*% power_estimator_S_Sigma2
  }
  
  ## Fast generation of Hankel matrices
  id <- 1:(m+1) + rep(0:m, each = m + 1)
  
  M <- matrix(s2[id], nrow = m + 1)
  
  # invM_solve = solve(M)
  
  invM_recursive = compute_M_inverse(m = m, all_tr0 = 1 / s2[1],
                                     all_tr = s2[-1], verbose = verbose - 2,
                                     mpfr = mpfr, precBits = precBits)
  
  alpha = invM_recursive %*% hm
  
  result = list(
    hm = hm,
    M = M,
    invM_recursive = invM_recursive,
    alpha = alpha
  )
  
  return (result)
}



loss_L2_ridge_oracle_higher_order_optimal <- function (
    S, t, Ip, p, m, Sigma, method_invM = "recursive", verbose, mpfr, precBits)
{
  if (verbose > 0){
    cat("t = ", t)
    if (verbose > 1){
      cat("\n")
    } else {
      cat(", ")
    }
  }
  
  iS_ridge <- solve(S + t * Ip)
  
  loss = tryCatch({
    estimatedM = compute_M_oracle(
      estimator_S = iS_ridge, Sigma = Sigma, m = m, Ip = Ip, p = p,
      verbose = verbose - 1, mpfr = mpfr, precBits = precBits)
    
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



loss_L2_MPR_oracle_higher_order_optimal <- function (
    S, t, Ip, p, m, Sigma, method_invM = "recursive", verbose, mpfr, precBits)
{
  if (verbose > 0){
    cat("t = ", t)
    if (verbose > 1){
      cat("\n")
    } else {
      cat(", ")
    }
  }
  
  iS_ridge <- solve(S + t * Ip)
  
  MPR_estimator <- iS_ridge - t * iS_ridge %*% iS_ridge
  
  loss = tryCatch({
    estimatedM = compute_M_oracle(
      estimator_S = MPR_estimator, Sigma = Sigma, m = m, Ip = Ip, p = p,
      verbose = verbose - 1, mpfr = mpfr, precBits = precBits)
    
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

