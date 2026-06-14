

#' 
#' @examples
#' 
#' for (m in 1:10){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = TRUE)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "Moore-Penrose", centeredCov = TRUE,
#'       method_invM = "recursive", verbose = 0)
#'       
#'   #print(precision_higher_order_shrinkage_Cent$M)
#'   cat("oracle$M = \n")
#'   print(oracle$M)
#'   
#'   cat("oracle$invM = \n")
#'   #print(oracle$invM_solve)
#'   print(oracle$invM_recursive)
#'   
#'   
#'   cat("estimated inv M = \n")
#'   print(precision_higher_order_shrinkage_Cent$invM)
#' }
#' 
#' for (m in 1:10){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = TRUE)
#'   
#'   oracle = oracle_higher_order_shrinkage(
#'       X = X, m = m, Sigma = Sigma, nameEstimator = "Moore-Penrose", centeredCov = TRUE,
#'       method_invM = "recursive", verbose = 0)
#'       
#'   #print(precision_higher_order_shrinkage_Cent$M)
#'   cat("oracle$hm = \n")
#'   print(oracle$hm)
#'   
#'   cat("estimated hm = \n")
#'   print(precision_higher_order_shrinkage_Cent$hm)
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
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- Moore_Penrose(X = X, centeredCov = centeredCov)
  iS_MP_ <- as.matrix(iS_MP)
  
  
  resultM <- compute_M_oracle(
    estimator_S = iS_MP_, Sigma = Sigma, m = m, Ip = Ip, p = p,
    verbose = verbose)
  
  result = list(
    M = resultM$M,
    hm = resultM$hm,
    # invM_solve = resultM$invM_solve,
    invM_recursive = resultM$invM_recursive
  )
  
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
  
  return (list(hm = hm, M = M,
               # invM_solve = invM_solve,
               invM_recursive = invM_recursive))
}

