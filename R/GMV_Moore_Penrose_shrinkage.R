
#' First-order shrinkage of the Moore-Penrose GMV portfolio towards a general target 
#' portfolio \eqn{\mathbf{b}}
#'
#' This function computes shrinked optimal portfolio weights given by
#' \deqn{\hat{\alpha}^*\mathbf{w}_{MP} + (1 - \hat{\alpha}^*) \mathbf{b}}
#' where \eqn{\hat{\alpha}^*} is given by
#' \deqn{
#' \hat{\alpha}^*=
#' \dfrac{\mathbf{b}^\top\mathbf{S}_n\mathbf{b}-
#' \frac{\hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma} \right)}
#' {\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}}
#' {\mathbf{b}^\top\mathbf{S}_n\mathbf{b} - 
#' 2\frac{\hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma} \right)}
#' {\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
#' +\frac{\hat{d}_3\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
#' {\hat{d}_1^2\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}},
#' } where
#' \deqn{
#' \hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma}\right) 
#' = \frac{1}{\hat{v}(0)}
#' \left[\frac{1}{\hat{v}(0)}\left(1-\hat{d}_0(0,\mathbf{1}\mathbf{b}^\top) \right) 
#' - \hat{d}_1(\mathbf{1}\mathbf{b}^\top) \right],
#'}
#' with \eqn{\hat{v}(0)}, \eqn{\hat{d}_1\left( \mathbf{1}\mathbf{b}^\top\right)}, 
#' \eqn{\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
#' \eqn{d_3(\frac{\mathbf{1}\mathbf{1}^\top}{p})}, and
#' \eqn{\hat{d}_0\left(0, \mathbf{1}\mathbf{b}^\top\right)} given in (8) and 
#' in (S.41), (S.44), and (S.45)  from the supplement of Bodnar and Parolya (2026).
#' The vector \eqn{\mathbf{w}_{MP}} are the optimal portfolio weights 
#' estimated as the plug-in of the Moore-Penrose estimate of the precision
#' matrix and \eqn{\mathbf{1}} is a vector of ones of size \eqn{p}.
#'
#'
# In the particular case of the shrinkage to the equally weighted portfolio,
# this function computes
# \deqn{\hat{\alpha}^*\times \mathbf{w}_{MP} + (1 - \hat{\alpha}^*)\
# \times \mathbf{1}/p}
# where \eqn{\hat{\alpha}^*} is given by
# \deqn{
# \hat{\alpha}^*=
# \frac{\frac{1}{p}\mathbf{1}^\top\mathbf{S}_n\mathbf{1}-
# \frac{\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}
# \boldsymbol{\Sigma} \right)}
# {\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}}
# {\frac{1}{p}\mathbf{1}^\top\mathbf{S}_n\mathbf{1} - 
# 2\frac{\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}
# \boldsymbol{\Sigma} \right)}
# {\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
# +\frac{\hat{d}_3\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
# {\hat{d}_1^2\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}},
# } where
# \deqn{
# \hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\boldsymbol{\Sigma}\right)
# = \frac{1}{\hat{v}(0)}
# \left[\frac{1}{\hat{v}(0)}\left(1-\hat{d}_0(0,
# \frac{\mathbf{1}\mathbf{1}^\top}{p}) \right) 
# - \hat{d}_1(\frac{\mathbf{1}\mathbf{1}^\top}{p}) \right],
#}
# with \eqn{\hat{v}(0)}, \eqn{\hat{d}_1\left(
# \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
# \eqn{\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
# \eqn{d_3(\frac{\mathbf{1}\mathbf{1}^\top}{p})}, and
# \eqn{\hat{d}_0\left(0, \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)} given
# in (8) and in (S.41), (S.44), and (S.45)  from the supplement of
# Bodnar and Parolya (2026).
# The vector \eqn{w_{MP}} are the optimal portfolio weights estimated as the
# plug-in of the Moore-Penrose estimate of the precision matrix and
# \eqn{\mathbf{1}} is a vector of ones of size \eqn{p}.
#
#'
#'
#' @template param-X
#' 
#' @param b shrinkage target. By default, the equally-weighted portfolio is used
#' as a target.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns An object of class `EstimatedPortfolioWeights`. The estimated
#' weights can be extracted as a numeric vector with `as.numeric()`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.1214/25-AOS2602}. ArXiv: \doi{10.48550/arXiv.2403.15792}.
#' 
#' @examples
#' set.seed(1)
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
#' # Compute GMV portfolio based on the Moore-Penrose inverse (no shrinkage)
#' GMV_MP = GMV_Moore_Penrose(X)
#' 
#' Losses(GMV_MP, Sigma)
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse with shrinkage
#' # towards the equally weighted portfolio
#' GMV_MP_shrink_eq = GMV_Moore_Penrose_shrinkage(X)
#' 
#' Losses(GMV_MP_shrink_eq, Sigma)
#' 
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse with shrinkage
#' # towards the true GMV portfolio
#' GMV_true = GMV_PlugIn(solve(Sigma))
#' 
#' GMV_MP_shrink_oracle = GMV_Moore_Penrose_shrinkage(X, b = GMV_true)
#' 
#' Losses(GMV_MP_shrink_oracle, Sigma)
#' 
#' 
#' @export
GMV_Moore_Penrose_shrinkage <- function(X, centeredCov = TRUE, b = NULL,
                                        verbose = 0){
  call_ = match.call()
  if (is.null(b)){
    if (verbose > 0){
      cat("Default target: equally weighted portfolio\n")
    }
    result = GMV_Moore_Penrose_shrinkage_eq(
      X = X, centeredCov = centeredCov, verbose = verbose, call_ = call_)
  } else {
    if (verbose > 0){
      cat("User-provided target portfolio\n")
    }
    result = GMV_Moore_Penrose_shrinkage_general(
      X = X, centeredCov = centeredCov, b = b, verbose = verbose, call_ = call_)
  }
  
  return (result)
}



GMV_Moore_Penrose_shrinkage_general <- function(X, centeredCov = TRUE, b = NULL,
                                                verbose = 0, call_ = NULL){
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Vector of ones of size p
  ones = rep(1, length = p)
  
  b = prepare_and_check_b(b = b, p = p)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(X = X, centeredCov = centeredCov))
  
  MP_portfolio = GMV_PlugIn(precisionMatrix = iS_MP)
  w_MP = as.numeric(MP_portfolio)
  
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
  
  
  result = new_estimated_portfolio_weights(
    estimate = w_ShMP,
    n = n,
    p = p,
    alpha_optimal = alp_ShMP,
    target = b,
    centeredCov = centeredCov,
    method = "Moore-Penrose shrinkage",
    call = call_
  )
  
  return (result)
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



GMV_Moore_Penrose_shrinkage_eq <- function(X, centeredCov = TRUE, verbose = 0,
                                           call_ = NULL){
  
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(X = X, centeredCov = centeredCov))
  
  MP_portfolio = GMV_PlugIn(precisionMatrix = iS_MP)
  w_MP = as.numeric(MP_portfolio)
  
  trS1 <- sum(diag(iS_MP)) / p
  trS2 <- sum(diag(iS_MP %*% iS_MP))/p
  trS3 <- sum(diag(iS_MP %*% iS_MP %*% iS_MP))/p
  trS4 <- sum(diag(iS_MP %*% iS_MP %*% iS_MP %*% iS_MP))/p
  
  bip<-matrix(rep(1,p),p,1)
  tbip<-t(bip)
  
  # No division by p here below
  
  bipiSbip  <- sum(tbip %*% iS_MP %*% bip)
  bipiS2bip <- sum(tbip %*% iS_MP %*% iS_MP %*% bip)
  bipiS3bip <- sum(tbip %*% iS_MP %*% iS_MP %*% iS_MP %*% bip)
  
  bipSbip   <- sum(tbip %*% S %*% bip)
  
  hv0 <- c_n * trS1
  
  d0   <- 1 - tbip %*% S %*% iS_MP %*% bip / p
  d1   <- bipiSbip / (p * c_n * trS2)
  d1_bSigma <- ( (1 - d0) / hv0 - d1) / hv0
  
  d3_first_term = bipiS3bip / (p * c_n^3 * trS2^3)
  d3_second_term = 2 * bipiSbip * trS3^2 / (c_n^3 * p * trS2^5)
  d3_third_term = (bipiS2bip + trS4 * bipiSbip) / (p * c_n^3 * trS2^4)
  
  d3 = d3_first_term + d3_second_term - d3_third_term
  
  if (verbose > 1){
    cat("  * d3_first_term  = ", d3_first_term, "\n")
    cat("  * d3_second_term = ", d3_second_term, "\n")
    cat("  * d3_third_term  = ", d3_third_term, "\n")
  }
  
  if (verbose > 1){
    cat("* hv0 = ", hv0, "\n")
    cat("* d0 = ", d0, "\n")
    cat("* d1 = ", d1, "\n")
    cat("* d1_bSigma = ", d1_bSigma, "\n")
    cat("* d3 = ", d3, "\n")
  }
  
  num = sum(bipSbip / p - d1_bSigma  / d1)
  den = sum(bipSbip / p - 2 * d1_bSigma / d1 + d3 / d1^2)
  
  alp_ShMP <- num / den
  
  if (verbose > 0){
    cat("num = ", num, "\n")
    cat("den = ", den, "\n")
    cat("alp_ShMP = ", alp_ShMP, "\n")
  }
  
  w_ShMP <- alp_ShMP * w_MP + (1 - alp_ShMP) * rep(1,p) / p
  
  result = new_estimated_portfolio_weights(
    estimate = w_ShMP,
    n = n,
    p = p,
    alpha_optimal = alp_ShMP,
    target = "equally weighted",
    centeredCov = centeredCov,
    method = "Moore-Penrose shrinkage",
    call = call_
  )
  
  return (result)
}

