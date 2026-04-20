

#' First-order shrinkage of the Ridge portfolio towards a general target 
#' portfolio \eqn{\mathbf{b}}
#'
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param b shrinkage target. By default, the equally-weighted portfolio is used
#' as a target.
#' 
#' @inheritParams cov_with_centering
#' 
#' @return a vector of size \eqn{p} of (estimated) optimal portfolio weights,
#' where \eqn{p} is the number of assets.
#' 
#' @references
#' Bodnar, T., Parolya, N., & Thorsén, E. (2024).
#' Two is better than one: Regularized shrinkage of large minimum variance
#' portfolios.
#' Journal of Machine Learning Research, 25(173), 1-32.
#' 
#' @examples
#' set.seed(1)
#' n = 50
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
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse (no shrinkage)
#' GMV_MP = GMV_Moore_Penrose(X)
#' 
#' Loss_GMV_MP = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_MP, Sigma = Sigma)
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse with shrinkage
#' # towards the equally weighted portfolio
#' GMV_MP_shrink_eq = GMV_Moore_Penrose_shrinkage(X)
#' 
#' Loss_GMV_MP_shrink_eq = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_MP_shrink_eq, Sigma = Sigma)
#' 
#' 
#' # Compute GMV portfolio based on the Ridge with shrinkage towards the equally
#' # weighted portfolio
#' GMV_ridge_shrink_eq = GMV_ridge_shrinkage(X)
#' 
#' Loss_GMV_ridge_shrink_eq = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_ridge_shrink_eq, Sigma = Sigma)
#'
#' cat("Loss of GMV_MP = ", Loss_GMV_MP, "\n")
#' cat("Loss of GMV_MP_shrink = ", Loss_GMV_MP_shrink_eq, "\n")
#' cat("Loss of GMV_ridge_shrink = ", Loss_GMV_ridge_shrink_eq, "\n")
#' 
#' @export
#' 
GMV_ridge_shrinkage <- function(X, centeredCov = TRUE, b = NULL,
                                verbose = 0, 
                                eps = 1/(10^6), upp = pi/2 - eps,
                                initialValue = 1.5){
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  
  # Vector of ones of size p
  ones = rep(1, length = p)
  
  if (is.null(b)){
    b = rep(1/p, length = p)
  } else if (length(b) != p){
    stop("'b' should be a vector of length 'p'.")
  }
  if (abs(sum(b) - 1) > 0.001){
    stop("The weights (b) should sum up to 1.")
  }
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Wrapper around the loss to optimize
  hL2R <- function(u){
    eta = tan(u)
    
    ridge_ = ridge(X = X, centeredCov = centeredCov, t = eta, verbose = verbose - 2,
                   method_inversion = "auto")
    
    iS_ridge <- as.matrix(ridge_)
    
    loss = loss_GMV_ridge_shrinkage(
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
  eta_opt <- tan(u_R)
  if (verbose > 0){
    cat("*  optimal eta =", eta_opt ,"\n")
  }
  
  ridge_ = ridge(X = X, centeredCov = centeredCov, t = eta_opt,
                 verbose = verbose - 2, method_inversion = "auto")
  
  alpha_opt = BNP2024_alpha_optimal(
    S = S, iS_ridge = as.matrix(ridge_), eta = eta_opt, ones = ones, b = b,
    p = p, c_n = c_n, verbose = verbose - 1)
  
  ridge_portfolio = GMV_PlugIn(ridge_)
  
  shrinked_portfolio = alpha_opt * ridge_portfolio + (1 - alpha_opt) * b
  
  return (shrinked_portfolio)
}



loss_GMV_ridge_shrinkage  <- function(S, iS_ridge, eta, ones, b, p, c_n, verbose){
  tr_iS_ridge = tr(iS_ridge)
  tr_iS_ridge2 = tr(iS_ridge %*% iS_ridge)
  vhat = BNP2024_estimator_vhat(tr_iS_ridge = tr_iS_ridge, p = p,
                                eta = eta, c_n = c_n)
  vhat_1_1 = BNP2024_estimator_vhat_1_1(
    tr_iS_ridge = tr_iS_ridge, tr_iS_ridge2 = tr_iS_ridge2,
    p = p, eta = eta, c_n = c_n, vhat = vhat)
  vhat_1_2 = BNP2024_estimator_vhat_1_2(vhat = vhat, vhat_1_1 = vhat_1_1,
                                        eta = eta)
  d1 = BNP2024_estimator_d1(iS_ridge = iS_ridge, b = b, ones = ones, eta = eta,
                            vhat = vhat)
  d2 = BNP2024_estimator_d2(iS_ridge = iS_ridge, eta = eta, vhat = vhat)
  
  # Term b^T S_n b
  bSb = t(b) %*% S %*% b
  numerator_term2 = d1 / ( (1 + eta) * sum(iS_ridge) )
  numerator = (bSb - numerator_term2) * (1 - numerator_term2 / bSb)
  
  denominator_term2 = 2 * numerator_term2
  denominator_term3_num = (1 - vhat_1_2) * d2
  denominator_term3_den = (sum(iS_ridge))^2 * (1 + eta)^2
  denominator_term3 = denominator_term3_num / denominator_term3_den
  denominator = bSb - denominator_term2 + denominator_term3
  
  # Note that in the paper, what is written is the opposite of the loss.
  # Therefore, in order to get the loss, we put a minus sign here below.
  result = - numerator / denominator
  
  if (verbose > 0){
    cat("loss = ", result, "\n")
    cat("* numerator = ", numerator, "\n")
    cat("  * b^T S_n b = ", bSb, "\n")
    cat("  * term2 = ", numerator_term2, "\n")
    cat("* denominator = ", denominator, "\n")
    cat("  * term2 = ", denominator_term2, "\n")
    cat("  * term3 = ", denominator_term3, "\n")
    cat("Other variables: \n")
    cat("* tr_iS_ridge = ", tr_iS_ridge, "\n")
    cat("* tr_iS_ridge2 = ", tr_iS_ridge2, "\n")
    cat("* vhat = ", vhat, "\n")
    cat("* vhat_1_1 = ", vhat_1_1, "\n")
    cat("* vhat_1_2 = ", vhat_1_2, "\n")
    cat("* d1 = ", d1, "\n")
    cat("* d2 = ", d2, "\n")
  }
  
  return (result)
}


# Theorem 2 page 8 of (Bodnar, T., Parolya, N., & Thorsén, E., 2024)
BNP2024_estimator_vhat <- function(tr_iS_ridge, p, eta, c_n){
  term_inside = 1 - eta * tr_iS_ridge / p
  result = 1 - c_n * tr_iS_ridge
  
  return (result)
}

# Theorem 2 page 8 of (Bodnar, T., Parolya, N., & Thorsén, E., 2024)
BNP2024_estimator_vhat_1_1 <- function(tr_iS_ridge, tr_iS_ridge2, p, eta,
                                       c_n, vhat){
  term_inside = tr_iS_ridge / p - eta * tr_iS_ridge2 / p
  
  result = vhat * c_n * term_inside
  
  return (result)
}

# Theorem 2 page 8 of (Bodnar, T., Parolya, N., & Thorsén, E., 2024)
BNP2024_estimator_vhat_1_2 <- function(vhat, vhat_1_1, eta){
  term_2 = 1 / vhat
  term_3 = eta * vhat_1_1 / vhat^2
  
  result = 1 - term_2 + term_3
  
  return (result)
}

# Inspired from formula (31) page 9 of
# (Bodnar, T., Parolya, N., & Thorsén, E., 2024)
BNP2024_estimator_d1 <- function(iS_ridge, b, ones, eta, vhat){
  matrix_product = as.numeric(t(b) %*% iS_ridge %*% ones)
  stopifnot(length(matrix_product) == 1)
  
  term_inside = 1 - eta * matrix_product
  result = ((1 + eta) / vhat) * term_inside
  
  return (result)
}


# Inspired from formula (32) page 9 of
# (Bodnar, T., Parolya, N., & Thorsén, E., 2024)
BNP2024_estimator_d2 <- function(iS_ridge, eta, vhat){
  sum1 = sum(iS_ridge)
  sum2 = sum(iS_ridge^2)
  
  term_inside = sum1 - eta * sum2
  result = ((1 + eta)^2 / vhat) * term_inside
  
  return (result)
}


BNP2024_alpha_optimal <- function(S, iS_ridge, eta, ones, b, p, c_n, verbose){
  tr_iS_ridge = tr(iS_ridge)
  tr_iS_ridge2 = tr(iS_ridge %*% iS_ridge)
  vhat = BNP2024_estimator_vhat(tr_iS_ridge = tr_iS_ridge, p = p,
                                eta = eta, c_n = c_n)
  vhat_1_1 = BNP2024_estimator_vhat_1_1(
    tr_iS_ridge = tr_iS_ridge, tr_iS_ridge2 = tr_iS_ridge2,
    p = p, eta = eta, c_n = c_n, vhat = vhat)
  vhat_1_2 = BNP2024_estimator_vhat_1_2(vhat = vhat, vhat_1_1 = vhat_1_1,
                                        eta = eta)
  d1 = BNP2024_estimator_d1(iS_ridge = iS_ridge, b = b, ones = ones, eta = eta,
                            vhat = vhat)
  d2 = BNP2024_estimator_d2(iS_ridge = iS_ridge, eta = eta, vhat = vhat)
  
  # Term b^T S_n b
  bSb = t(b) %*% S %*% b
  numerator_term2 = d1 / ( (1 + eta) * sum(iS_ridge) )
  numerator = bSb - numerator_term2
  
  denominator_term2 = 2 * numerator_term2
  denominator_term3_num = (1 - vhat_1_2) * d2
  denominator_term3_den = (sum(iS_ridge))^2 * (1 + eta)^2
  denominator_term3 = denominator_term3_num / denominator_term3_den
  denominator = bSb - denominator_term2 + denominator_term3
  
  result = numerator / denominator
  
  if (verbose > 0){
    cat("alpha_optimal = ", result, "\n")
    cat("* numerator = ", numerator, "\n")
    cat("  * b^T S_n b = ", bSb, "\n")
    cat("  * term2 = ", numerator_term2, "\n")
    cat("* denominator = ", denominator, "\n")
    cat("  * term2 = ", denominator_term2, "\n")
    cat("  * term3 = ", denominator_term3, "\n")
    cat("Other variables: \n")
    cat("* tr_iS_ridge = ", tr_iS_ridge, "\n")
    cat("* tr_iS_ridge2 = ", tr_iS_ridge2, "\n")
    cat("* vhat = ", vhat, "\n")
    cat("* vhat_1_1 = ", vhat_1_1, "\n")
    cat("* vhat_1_2 = ", vhat_1_2, "\n")
    cat("* d1 = ", d1, "\n")
    cat("* d2 = ", d2, "\n")
  }
  
  return (result)
}
  