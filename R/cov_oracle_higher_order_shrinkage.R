
#' Oracle higher order shrinkage for estimation of the covariance matrix
#'
#' This function computes the oracle estimator of the covariance
#' matrix using the polynomial
#' \deqn{
#' \mathbf{S}_{n;HOS}^{or}(m)
#' =\alpha_0(m)\mathbf{I}_p
#' +\sum_{j=1}^m\alpha_j(m)\mathbf{S}_n^j,
#' }
#' where \eqn{\boldsymbol{\alpha}(m)
#' =(\alpha_0(m),\alpha_1(m),\ldots,\alpha_m(m))^\top}
#' is given by
#' \deqn{
#' \boldsymbol{\alpha}(m)
#' =\mathbf{M}_n(m)^{-1}\mathbf{m}_n(m)
#' }
#' with
#' \deqn{
#' \mathbf{m}_n(m)=
#' \begin{pmatrix}
#' \frac{1}{p}\operatorname{tr}[\boldsymbol{\Sigma}]\\
#' \frac{1}{p}\operatorname{tr}[\mathbf{S}_n\boldsymbol{\Sigma}]\\
#' \vdots\\
#' \frac{1}{p}\operatorname{tr}[\mathbf{S}_n^m\boldsymbol{\Sigma}]
#' \end{pmatrix}
#' \quad \text{and}\quad
#' \mathbf{M}_n(m)=
#' \begin{pmatrix}
#' 1 & \hat{t}_1 & \ldots & \hat{t}_m\\
#' \hat{t}_1 & \hat{t}_2 & \ldots & \hat{t}_{m+1}\\
#' \vdots & \vdots & \ddots & \vdots\\
#' \hat{t}_m & \hat{t}_{m+1} & \ldots & \hat{t}_{2m}
#' \end{pmatrix},
#' }
#' where \eqn{\hat{t}_j=p^{-1}\operatorname{tr}[\mathbf{S}_n^j]},
#' \eqn{\mathbf{S}_n} is the sample covariance matrix,
#' \eqn{\boldsymbol{\Sigma}} is the population covariance matrix,
#' and \eqn{\mathbf{I}_p} is the identity matrix of size \eqn{p}.
#' The inverse of \eqn{\mathbf{M}_n(m)} is understood whenever
#' this matrix is nonsingular. This procedure ensures that the loss 
#'\eqn{||\mathbf{S}_{n;HOS}^{or}(m)-\boldsymbol{\Sigma}||^2_F}
#' is minimized with probability one for a finite sample.
#' 
#' @template param-X
#' 
#' @inheritParams cov_with_centering
#' @param Sigma true covariance matrix
#' 
#' @template param-m
#' @template param-mpfr
#' 
#' @template returns-EstimatedCovarianceMatrix
#' 
#' @examples
#' p = 300
#' n = 600
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' 
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' estimatedCov_sample = cov(X)
#' estimatedCov_analytical = cov_analytical_NL_shrinkage(X)
#' estimatedCov_QIS = cov_quadratic_inverse_shrinkage(X)
#' 
#' # We now compare the distances between the truth and the estimators.
#' Losses(estimatedCov_sample, Sigma, type = "covariance")
#' Losses(estimatedCov_analytical, Sigma)
#' Losses(estimatedCov_QIS, Sigma)
#' 
#' for (m in 0:5){
#'   estimatedCov_shrink_higher = cov_higher_order_shrinkage(X, m = m)
#'   
#'   loss = LossFrobenius2(estimatedCov_shrink_higher, Sigma)
#'   cat("m = ", m, ", loss = ", loss, "\n")
#'   
#'   print(coef(estimatedCov_shrink_higher))
#'   
#'   estimatedCov_shrink_higher_oracle = 
#'     cov_oracle_higher_order_shrinkage(X, m = m, Sigma = Sigma)
#'   
#'   loss = LossFrobenius2(estimatedCov_shrink_higher_oracle, Sigma)
#'   cat("m = ", m, ", loss = ", loss, "\n")
#'   
#'   print(coef(estimatedCov_shrink_higher_oracle))
#'   cat("\n")
#' }
#' 
#' oracle_Estimator = cov_NL_oracle(X, Sigma = Sigma)
#' Losses(oracle_Estimator, Sigma)
#' 
#' 
#' @export
cov_oracle_higher_order_shrinkage  <- function(
    X, centeredCov = TRUE, m, Sigma, verbose = 0,
    mpfr = FALSE, precBits = 2^16, warn_large_m = TRUE)
{
  call_ = match.call()
  
  if (verbose > 0){
    cat("Starting `cov_oracle_higher_order_shrinkage`...\n")
  }
  check_Rmpfr(mpfr)
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  check_large_m(warn = warn_large_m, p = p, m = m)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  if (verbose > 0){
    cat("*  m = ", m, "\n")
  }
  
  list_power_S = list()
  power_S = Ip
  
  for (k in seq_len(2 * m)){
    if (k == 1){
      power_S = S
    } else {
      power_S = power_S %*% S
    }
    list_power_S[[k]] <- power_S
  }
  
  estimatedM = compute_M_covariance_oracle(
    m = m, c_n = c_n, p = p, S = S, Sigma = Sigma, verbose = verbose,
    list_power_S = list_power_S, mpfr = mpfr, precBits = precBits)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = estimatedM$invM %*% estimatedM$hm
  if (verbose > 0){
    cat("Optimal alpha: \n")
    print(alpha)
  }
  
  estimated_covariance_matrix = alpha[1] * Ip
  power_S = Ip
  
  # Here the case m = 0 is handled since seq_len(0)=numeric(0)
  # and therefore skips the loop.
  for (k in seq_len(m)){
    estimated_covariance_matrix = estimated_covariance_matrix + 
      alpha[k + 1] * list_power_S[[k]]
  }
  
  result = new_estimated_covariance_matrix(
    estimate = estimated_covariance_matrix,
    n = n,
    p = p,
    centeredCov = centeredCov,
    method = "cov_oracle_higher_order_shrinkage",
    call = call_,
    invM = estimatedM$invM,
    hm = estimatedM$hm,
    alpha = alpha
  )
  
  return (result)
}


compute_M_covariance_oracle <- function(
    m, c_n, p, S, Sigma, verbose, list_power_S, mpfr, precBits)
{
  hm = rep(NA, m+1)
  
  for (k in seq_len(m + 1)){
    if (k == 1){
      power_Sigma_S = Sigma
    } else {
      power_Sigma_S = S %*% power_Sigma_S
    }
    
    hm[k] = tr(power_Sigma_S) / p
  }
  
  if (m == 0){
    
    result = list(invM = 1, hm = hm)
    return (result)
  }
  
  # Now m is at least 1
  
  all_tr = compute_trace_power_matrix(list_power_S = list_power_S,
                                      m = 2 * m, verbose = verbose - 1)
  
  invM = compute_M_inverse(
    m = m,
    all_tr0 = 1, # For covariance matrices, the initialization is at 1.
    all_tr = all_tr / p, verbose = verbose,
    mpfr = mpfr, precBits = precBits)
  
  return (list(invM = invM, hm = hm))
}

