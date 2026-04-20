

#' 
#' @export
GMV_MPR_shrinkage <- function(X, centeredCov = TRUE, b = NULL,
                                verbose = 0, 
                                eps = 1/(10^6), upp = pi/2 - eps,
                                initialValue = 1.5){
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  
  # Vector of ones of size p
  ones = rep(1, length = p)
  
  b = prepare_and_check_b(b = b, p = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  
  # Wrapper around the loss to optimize
  hL2R <- function(u){
    eta = tan(u)
    
    # ridge_ = ridge(X = X, centeredCov = centeredCov, t = eta, verbose = verbose - 2,
    #                method_inversion = "auto")
    # 
    # iS_ridge <- as.matrix(ridge_)
    
    loss = loss_GMV_MPR_shrinkage(
      S = S, iS_ridge = iS_ridge, eta = eta, ones = ones, b = b, p = p,
      c_n = c_n, verbose = verbose - 3)
    
    if (verbose > 1){
      cat("eta = ", eta, "; loss = ", loss, ". ")
    }
    
    return(loss)
  }
  
  control = list(trace = if(verbose > 2){6} else {0},
                 factr = 1e8
                 # ndeps = 0.01
  )
  
  hL2R_max <- stats::optim(par = initialValue, fn = hL2R,
                           lower = eps, upper = upp,
                           method = "L-BFGS-B", control = control)
  
  u_R <- hL2R_max$par
  t_opt <- tan(u_R)
  if (verbose > 0){
    cat("*  optimal t =", t_opt ,"\n")
  }
  
  MPR_ = MPR(X = X, centeredCov = TRUE, t = t_opt, verbose = verbose - 3)
  
  alpha_opt = alpha_optimal_GMV_MPR()
  
  MPR_portfolio = GMV_PlugIn(MPR_)
  
  shrinked_portfolio = alpha_opt * MPR_portfolio + (1 - alpha_opt) * b
  
  return (shrinked_portfolio)
}


loss_GMV_ridge_shrinkage  <- function(S, iS_ridge, t, ones, b, p, c_n, verbose){

  
}


alpha_optimal_GMV_MPR <- function(){
  
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(X = X, centeredCov = centeredCov))
  
  
  trS1 <- sum(diag(iS_MP)) / p
  trS2 <- sum(diag(iS_MP %*% iS_MP))/p
  trS3 <- sum(diag(iS_MP %*% iS_MP %*% iS_MP))/p
  trS4 <- sum(diag(iS_MP %*% iS_MP %*% iS_MP %*% iS_MP))/p
  
  ones = matrix(data = 1, nrow = p, ncol = 1)
  bip  = matrix(data = b, nrow = p, ncol = 1)
  tones <- t(ones)
  tbip  <- t(bip)
  
  ones_iS_ones  <- sum(tones %*% iS_MP %*% ones)
  ones_iS2_ones <- sum(tones %*% iS_MP %*% iS_MP %*% ones)
  ones_iS3_ones <- sum(tones %*% iS_MP %*% iS_MP %*% iS_MP %*% ones)
  
  bipSbip <- sum(tbip%*%S%*%bip)
  
  hv0 <- c_n * trS1
  
  d0 <- 1 - tbip %*% S %*% iS_MP %*% ones
  
  d1   <- estimator_GMV_d1(iS = iS_MP,
                           Theta = matrix(1/p, nrow = p, ncol = p),
                           c_n = c_n, p = p, verbose = verbose - 1)
  
  d1_b <- estimator_GMV_d1(iS = iS_MP,
                           Theta = bip %*% tones,
                           c_n = c_n, p = p, verbose = verbose - 1)
  
  # Note that you cannot use the function `estimator_GMV_d1` because
  # Sigma is unknown, so this is not a true function. We must use the expression:
  d1_bSigma <- ((1 - d0) / hv0 - d1_b) / hv0
  
  d3 <- estimator_GMV_d3(p = p, iS = iS_MP,
                         Theta = matrix(1/p, nrow = p, ncol = p),
                         c_n = c_n,
                         verbose = verbose - 1)
  
  # d3 <- (ones_iS_ones / trS2^3 + 2 * ones_iS_ones * trS3^2 / p / trS2^5 - 
  #        (ones_iS2_ones + trS4 * ones_iS_ones) / trS2^4) / c_n^3
  # d3 <- d3 / p
  
  if (verbose > 1){
    cat("* hv0 = ", hv0, "\n")
    cat("* d0 = ", d0, "\n")
    cat("* d1 = ", d1, "\n")
    cat("* d1_b = ", d1_b, "\n")
    cat("* d1_bSigma = ", d1_bSigma, "\n")
    cat("* d3 = ", d3, "\n")
  }
  
  num = sum(p * bipSbip - d1_bSigma / d1)
  den = sum(p * bipSbip - 2 * d1_bSigma / d1 + d3/d1^2)
  
  alp_ShMP <- num / den
  
  if (verbose > 0){
    cat("num = ", num, "\n")
    cat("den = ", den, "\n")
    cat("alp_ShMP = ", alp_ShMP, "\n")
  }
  
  
}


