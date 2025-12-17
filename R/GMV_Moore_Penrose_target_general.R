

#' First-order shrinkage of the Moore-Penrose portfolio towards a general target portfolio
#'
#' This function computes
#' \deqn{\alpha \times w_{MP} + (1 - \alpha) \times \mathbf{b}}
#' where \eqn{\alpha} is a carefully chosen coefficient, \eqn{w_{MP}} is the
#' vector of optimal portfolio weights estimated as the plug-in of the Moore-Penrose
#' estimate of the precision matrix and \eqn{\mathbf{b}} is the target.
#'
#'
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param b target portfolio. This is a vector of size \code{p} which sums to 1.
#' By default this is the equally weighted portfolio.
#' 
#' @return a vector of size \eqn{p} of (estimated) optimal portfolio weights,
#' where \eqn{p} is the number of assets.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2024).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \link{https://doi.org/10.48550/arXiv.2403.15792}
#' 
#' 
#' @examples
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
#' GMV_MP_shrinkage_Cent = 
#'   GMV_Moore_Penrose_target_eq(Y = t(X), centeredCov = TRUE)
#'   
#' GMV_MP_shrinkage_Cent_new = 
#'   GMV_Moore_Penrose_target_general(Y = t(X), centeredCov = TRUE)
#'   
#' # This should coincide with GMV_MP_shrinkage_Cent since by default this is
#' # shrinked to the equally weighted portfolio.
#' 
#' outOfSampleVariance = t(GMV_MP_shrinkage_Cent) %*% Sigma %*% GMV_MP_shrinkage_Cent
#' 
#' ones = rep(1, length = p)
#' V_GMV = 1 / ( t(ones) %*% solve(Sigma) %*% ones)
#' 
#' Loss_GMV_Moore_Penrose_target_eq = (outOfSampleVariance - V_GMV) / V_GMV
#' 
#' GMV_MP_Cent = GMV_Moore_Penrose(Y = t(X), centeredCov = TRUE)
#' outOfSampleVariance = t(GMV_MP_Cent) %*% Sigma %*% GMV_MP_Cent
#' 
#' Loss_GMV_Moore_Penrose = (outOfSampleVariance - V_GMV) / V_GMV
#' 
#' # Shrinkage helps to reduce the loss
#' stopifnot(Loss_GMV_Moore_Penrose_target_eq < Loss_GMV_Moore_Penrose)
#' 
#' 
#' # We now compare with the true target
#' perfect_GMV_portfolio = GMV_PlugIn(solve(Sigma))
#' 
#' GMV_MP_shrinkage_Cent_oracle = 
#'   GMV_Moore_Penrose_target_general(Y = t(X), centeredCov = TRUE,
#'                                    b = perfect_GMV_portfolio)
#'                                    
#' outOfSampleVariance = 
#'   t(GMV_MP_shrinkage_Cent_oracle) %*% Sigma %*% GMV_MP_shrinkage_Cent_oracle
#' 
#' outOfSampleVarianceGMV = 
#'   t(perfect_GMV_portfolio) %*% Sigma %*% perfect_GMV_portfolio
#' 
#' ones = rep(1, length = p)
#' V_GMV = 1 / ( t(ones) %*% solve(Sigma) %*% ones)
#' 
#' Loss_GMV_MP_shrinkage_Cent_oracle = (outOfSampleVariance - V_GMV) / V_GMV
#' 
#' cat("GMV_Moore_Penrose_target_eq:", Loss_GMV_Moore_Penrose_target_eq, "\n")
#' cat("GMV_Moore_Penrose:", Loss_GMV_Moore_Penrose, "\n")
#' cat("GMV_Moore_Penrose_target_oracle:", Loss_GMV_MP_shrinkage_Cent_oracle, "\n")
#' 
#' 
#' @export
GMV_Moore_Penrose_target_general <- function(Y, centeredCov = TRUE, b = NULL,
                                             verbose = 2){
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
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

  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / (n-1)
    
    # We remove the last eigenvector because the eigenvalues are sorted
    # in decreasing order.
    Hn = eigen(Jn)$vectors[, -n]
    Ytilde = Y %*% Hn
    
    # Inverse companion covariance
    iYtilde <- solve(t(Ytilde) %*% Ytilde / (n-1))
    
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
  
  w_MP = GMV_PlugIn(estimatedPrecisionMatrix = iS_MP)
  
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
  
  
  w_ShMP <- alp_ShMP * w_MP + (1 - alp_ShMP) * b
  
  return (w_ShMP)
}


#' @param iS the Moore-Penrose inverse of the sample covariance matrix
#' 
#' @noRd
estimator_GMV_d3 <- function(p, iS, Theta, c_n, verbose){
  
  iS2 = iS %*% iS
  iS3 = iS2 %*% iS
  iS4 = iS3 %*% iS
  
  first_term_num = tr(iS3 %*% Theta)
  first_term_den = c_n^3 * (tr(iS2) / p)^3
  first_term = first_term_num / first_term_den
  
  second_term_num = 2 * (tr(iS3) / p)^2 * tr(iS %*% Theta)
  second_term_den = c_n^3 * (tr(iS2) / p)^5
  second_term = second_term_num / second_term_den
  
  third_term_num = tr(iS2 %*% Theta) + tr(iS4) * tr(iS %*% Theta) / p
  third_term_den = c_n^3 * (tr(iS2) / p)^4
  third_term = third_term_num / third_term_den
  
  if (verbose > 0){
    cat("  * d3_first_term  = ", first_term, "\n")
    cat("  * d3_second_term = ", second_term, "\n")
    cat("  * d3_third_term  = ", third_term, "\n")
  }
  
  result = first_term + second_term - third_term
  
  return (result)
}


estimator_GMV_d1 <- function(iS, Theta, c_n, p, verbose){
  
  num = tr(iS %*% Theta)
  den = c_n * (1/p) * tr(iS %*% iS)
  
  result = num / den
  return (result)
}
