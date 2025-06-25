

#' Ridge higher order shrinkage with optimal t
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' @param t penalization parameter
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
#' FrobeniusNorm2 <- function(M){sum(diag(M %*% t(M)))}
#' 
#' precision_MoorePenrose_Cent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = FALSE)
#'
#' FrobeniusLoss2 <- function(M){FrobeniusNorm2(M %*% Sigma - diag(p) ) / p}
#'
#' print(FrobeniusLoss2(precision_MoorePenrose_Cent))
#' print(FrobeniusLoss2(precision_MoorePenrose_NoCent))
#' 
#' for (m in 1:5){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       ridge_higher_order_shrinkage_optimal(Y = t(X), m = m, centeredCov = TRUE)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       ridge_higher_order_shrinkage_optimal(Y = t(X), m = m, centeredCov = FALSE)
#'       
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_Cent$estimated_precision_matrix))
#'   
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_NoCent$estimated_precision_matrix))
#' }
#' 
#' 
#' 
#' @export
#' 
ridge_higher_order_shrinkage_optimal <- function(Y, m, centeredCov, interval = c(0, 50))
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    c_n = p / (n-1)
    
  } else {
    S <- Y %*% t(Y)/n
    
    c_n = p / n
  }
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  
  estimatedLoss <- function (t){
    
    # Regularized sample covariance matrix (Tikhonov regularization)
    S_t <- S + t * diag(nrow = p)
    
    S_t_inverse <- solve(S_t)
    
    loss = tryCatch({
      estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                               S_t_inverse = S_t_inverse,
                               t = t)
      
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
  
  result_optim <- optimize(f = estimatedLoss, interval = interval)
  
  optimal_t = result_optim$minimum
  
  
  # ============================================================================
  # We now compute the estimator using this optimal_t that was found.
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  S_t <- S + optimal_t * diag(nrow = p)
  
  S_t_inverse <- solve(S_t)
  
  estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                           S_t_inverse = S_t_inverse,
                           t = optimal_t)
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  
  result = alpha[1] * Ip
  power_S_t_inverse = Ip
  
  for (k in 1:m){
    power_S_t_inverse = power_S_t_inverse %*% S_t_inverse
    result = result + alpha[k + 1] * power_S_t_inverse
  }
  
  return(list(
    estimated_precision_matrix = result,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v,
    t = optimal_t,
    estimated_loss = result_optim$objective
  ) )
}

