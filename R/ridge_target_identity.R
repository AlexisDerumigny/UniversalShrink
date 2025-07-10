


#' Ridge with target set to the identity
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @examples
#' 
#' n = 100
#' p = 5 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
#' precision_ridge_target_Cent = 
#'     ridge_target_identity_optimal(Y = t(X), centeredCov = TRUE)
#'     
#' precision_ridge_target_NoCent = 
#'     ridge_target_identity_optimal(Y = t(X), centeredCov = FALSE)
#' 
#' FrobeniusNorm2 <- function(M){sum(diag(M %*% t(M)))}
#' 
#' FrobeniusNorm2(precision_ridge_target_Cent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_ridge_target_NoCent %*% Sigma - diag(p) ) / p
#' 
#' 
#' @export
ridge_target_identity_optimal <- function (Y, centeredCov){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    # We remove the last eigenvector because the eigenvalues are sorted
    # in decreasing order.
    Hn = eigen(Jn)$vectors[, -n]
    Ytilde = Y %*% Hn
    
    # Inverse companion covariance
    iYtilde <- solve(t(Ytilde) %*% Ytilde / (n-1) )
    
    # Moore-Penrose inverse
    iS_MP <- Ytilde %*% iYtilde %*% iYtilde %*% t(Ytilde) / (n-1)
    
    c_n = p / (n-1)
  } else {
    S <- Y %*% t(Y)/n
    
    # Inverse companion covariance
    iY <- solve(t(Y) %*% Y / n)
    
    # Moore-Penrose inverse
    iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
    
    c_n = p / n
  }
  
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
  
  # TODO: provide this as an option for the user
  initialValue = 1.5
  eps <- 1/(10^6)
  upp <- pi/2 - eps
  
  hL2R_max <- optim(par = initialValue, fn = hL2R,
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


#' @rdname ridge_target_identity_optimal
#' @export
ridge_target_identity_semioptimal <- function (Y, centeredCov, t, verbose = 2){
  
  if (verbose){
    cat("Starting `ridge_target_identity_semioptimal`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  if (verbose){
    cat("*  n = ", n, "\n")
    cat("*  p = ", p, "\n")
    cat("*  t = ", t, "\n")
  }
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    if (verbose){
      cat("*  centered case\n")
    }
    
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    # We remove the last eigenvector because the eigenvalues are sorted
    # in decreasing order.
    Hn = eigen(Jn)$vectors[, -n]
    Ytilde = Y %*% Hn
    
    # Inverse companion covariance
    iYtilde <- solve(t(Ytilde) %*% Ytilde / (n-1) )
    
    # Moore-Penrose inverse
    iS_MP <- Ytilde %*% iYtilde %*% iYtilde %*% t(Ytilde) / (n-1)
    
    c_n = p / (n-1)
  } else {
    if (verbose){
      cat("*  non-centered case\n")
    }
    
    S <- Y %*% t(Y)/n
    
    # Inverse companion covariance
    iY <- solve(t(Y) %*% Y / n)
    
    # Moore-Penrose inverse
    iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
    
    c_n = p / n
  }
  
  if (verbose){
    cat("*  c_n = ", c_n, "\n\n")
  }
  
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
  
  
  if (verbose){
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
  
  M = c(q2                                  , d0Sig2_t1 / t,
        d0Sig2_t1 / t                       , d0Sig2_t1 / t^2 + hvprt1 * d1Sig2_t1 / t) |>
    matrix(ncol = 2, byrow = TRUE)
  
  if (verbose){
    cat("Term of M[2,2] computed as:\n")
    cat("*  d0Sig2_t1 / t^2: ", d0Sig2_t1 / t^2, "\n")
    cat("*  hvprt1 * d1Sig2_t1 / t: ", hvprt1 * d1Sig2_t1 / t, "\n\n")
  }
  
  det_M = den_ShRt1 / t
  
  inv_M = c(d0Sig2_t1 / t^2 + hvprt1 * d1Sig2_t1 / t, - d0Sig2_t1 / t,
            - d0Sig2_t1 / t                         , q2 ) |>
    matrix(ncol = 2, byrow = TRUE) / 
    det_M
  
  
  alpha <- num_a_ShRt1 / den_ShRt1
  beta <- num_b_ShRt1 / den_ShRt1
  
  iS_ShRt1 <- alpha * iS_ridge + beta * Ip
  
  if (verbose){
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


#' @rdname ridge_target_identity_optimal
#' @export
ridge_target_identity <- function (Y, centeredCov, t, alpha, beta){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
  } else {
    S <- Y %*% t(Y)/n
  }
  
  
  iS_ridge <- solve(S + t * Ip)
  
  iS_ShRt1 <- alpha * iS_ridge + beta * Ip
  
  result = list(
    estimated_precision_matrix = iS_ShRt1
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

