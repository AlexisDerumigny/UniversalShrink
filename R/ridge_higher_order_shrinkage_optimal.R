

#' 
#' 
#' @rdname ridge_higher_order_shrinkage
#' @export
#' 
ridge_higher_order_shrinkage_optimal <- function(Y, m, centeredCov = TRUE,
                                                 interval = c(0, 50),
                                                 verbose = 0)
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  c_n <- concentration_ratio(n = n, p = p, centeredCov = centeredCov,
                             verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  
  estimatedLoss <- function (t){
    
    # Regularized sample covariance matrix (Tikhonov regularization)
    S_t <- S + t * diag(nrow = p)
    
    S_t_inverse <- solve(S_t)
    
    loss = tryCatch({
      estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                               S_t_inverse = S_t_inverse,
                               t = t, verbose = 0)
      
      loss = 1 - t(estimatedM$hm) %*% solve(estimatedM$M) %*% estimatedM$hm
    }, error = function(e){e}
    )
    
    if (inherits(loss, "simpleError") || !is.finite(loss)){
      loss <- .Machine$double.xmax
    }
    
    return (loss)
  }
  
  # initialValue = 1
  # eps <- 1/(10^6)
  # upp <- pi/2 - eps
  
  result_optim <- stats::optimize(f = estimatedLoss, interval = interval)
  
  optimal_t = result_optim$minimum
  
  
  # ============================================================================
  # We now compute the estimator using this optimal_t that was found.
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  S_t <- S + optimal_t * diag(nrow = p)
  
  S_t_inverse <- solve(S_t)
  
  estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                           S_t_inverse = S_t_inverse,
                           t = optimal_t, verbose = verbose)
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  
  estimated_precision_matrix = alpha[1] * Ip
  power_S_t_inverse = Ip
  
  for (k in 1:m){
    power_S_t_inverse = power_S_t_inverse %*% S_t_inverse
    estimated_precision_matrix = estimated_precision_matrix + 
      alpha[k + 1] * power_S_t_inverse
  }
  
  result = list(
    estimated_precision_matrix = estimated_precision_matrix,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v,
    t = optimal_t,
    estimated_loss = result_optim$objective
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

