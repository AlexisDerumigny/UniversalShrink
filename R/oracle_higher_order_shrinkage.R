

#' Oracle estimator of the precision matrix
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
#' for (m in 1:20){
#'   cat("\nm = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = Moore_Penrose_higher_order_shrinkage(
#'       X, m = m, centeredCov = TRUE, verbose = 0)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "Moore-Penrose",
#'       centeredCov = TRUE,method_invM = "recursive", verbose = 0)
#'       
#'   print(LossInverseFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   print(LossInverseFrobenius2(oracle, Sigma))
#' }
#' 
#' 
#' @export
oracle_higher_order_shrinkage <- function(
    X, m, Sigma, nameEstimator = "Moore-Penrose", centeredCov = TRUE,
    method_invM = "recursive", verbose = 0) {
  
  call_ = match.call()
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
    
  } else if (nameEstimator == "ridge") {
    
    
  } else if (nameEstimator == "MPR") {
    
  } else {
    stop("nameEstimator = '", nameEstimator, "' is not valid.",
         "It should be one of the following: 'Moore-Penrose', 'ridge', 'MPR'.")
  }
  
  
  resultM <- compute_M_oracle(
    estimator_S = estimated_precision_matrix, Sigma = Sigma,
    m = m, Ip = Ip, p = p, verbose = verbose)
  
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
    method = "Oracle higher-order shrinkage",
    nameBaselinEstimator = nameEstimator,
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

compute_M_oracle <- function(estimator_S, Sigma, m, Ip, p, verbose) {
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
                                     all_tr = s2[-1], verbose = verbose - 2)
  
  alpha = invM_recursive %*% hm
  
  result = list(
    hm = hm,
    M = M,
    invM_recursive = invM_recursive,
    alpha = alpha
  )
  
  return (result)
}

