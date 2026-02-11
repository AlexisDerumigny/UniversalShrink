
#' First-order shrinkage of the Moore-Penrose inverse towards a fixed target
#' 
#' Following Bodnar and Parolya (2006), the shrinkage estimator
#'  for the precision matrix using the Moore-Penrose inverse of the sample 
#'  covariance matrix \eqn{\mathbf{S}_n} is computed as
#' \deqn{
#' \widehat{\boldsymbol{\Pi}}_{MP}=\hat{\alpha}_{MP}^*\mathbf{S}_n^{+}+
#' \hat{\beta}_{MP}^*\boldsymbol{\Pi}_0\,,
#'} where \eqn{\mathbf{S}^+_n} denotes the Moore-Penrose inverse the sample
#' covariance matrix, \eqn{\boldsymbol{\Pi}_0} is the shrinkage target 
#' (\eqn{\boldsymbol{\Pi}_0=\mathbf{I}_p}, i.e., identity matrix by default) and
#'   \eqn{\hat{\alpha}_{MP}^*} and \eqn{\hat{\beta}_{MP}^*} are the optimal 
#'   shrinkage intensities (in sense of minimizing asymptotically 
#'   \eqn{||\widehat{\boldsymbol{\Pi}}_{MP}\boldsymbol{\Sigma}-\mathbf{I}_p||^2_F}) given by 
#' \deqn{
#' \hat{\alpha}_{MP}^* =
#'  \dfrac{
#'    \hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}\right)
#'    \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#'    -
#'      \hat{d}_1\left(0,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#'    \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#'  }{
#'    -\dfrac{1}{\hat{h}_2}
#'    \left(
#'      \hat{d}_2\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#'      -
#'        \hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#'      \dfrac{\hat{h}_3}{\hat{h}_2}
#'      \right)
#'    \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#'    -
#'      \dfrac{1}{\hat{h}_2}
#'    \hat{d}_1^2\left(0,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#'  }
#'} and
#' \deqn{
#' \hat{\beta}_{MP}^* =
#' \dfrac{
#' \left(
#' \hat{d}_2\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' -
#' \hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \dfrac{\hat{h}_3}{\hat{h}_2}
#' \right)
#' \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#' +
#' \hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{d}_1\left(0,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }{
#' \left(
#' \hat{d}_2\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' -
#' \hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \dfrac{\hat{h}_3}{\hat{h}_2}
#' \right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' +
#' \hat{d}_1^2\left(0,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }\,,
#' }
#' where \eqn{\hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}\right)}, 
#' \eqn{\hat{d}_1\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)}, 
#' \eqn{\hat{d}_1\left(0,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)},
#' \eqn{\hat{d}_2\left(\frac{1}{p}\boldsymbol{\Sigma}^2\right)} are defined in Theorem 3.2
#' in Bodnar and Parolya (2026), while the rest \eqn{\hat{v}(0)}, \eqn{\hat{h}_2}, \eqn{\hat{h}_3},
#' \eqn{\hat{d}_1(\frac{1}{p}\mathbf{I}_p)},
#' \eqn{\hat{d}_1(\frac{1}{p}\boldsymbol{\Pi}_0)},
#' \eqn{\hat{d}_2(\frac{1}{p}\mathbf{I}_p)},
#' \eqn{\hat{d}_0(0,\frac{1}{p}\boldsymbol{\Pi}_0)},
#' \eqn{\hat{q}_1(\frac{1}{p}\boldsymbol{\Pi}_0)}
#' and
#'\eqn{\hat{q}_2(\frac{1}{p}\boldsymbol{\Pi}_0^2)}
#' are given in the supplement to Bodnar and Parolya (2026).
#' 
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param Pi0 prior of the precision matrix. This a `p` by `p` matrix, used as
#' a target for the shrinkage. Default value is the identity matrix of size `p`.
#' As an advice, it should be a symmetric positive-definite matrix, but this is
#' not checked for.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix, i.e. the inverse of the covariance matrix).
#' 
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
#' 
#' @examples
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
#' precision_MoorePenrose_Cent = Moore_Penrose_target(X, centeredCov = TRUE)
#'    
#' precision_MoorePenrose_NoCent = Moore_Penrose_target(X, centeredCov = FALSE)
#' 
#' precisionTrue = solve(Sigma)
#' 
#' estimatedCov_NLshrink = cov_analytical_NL_shrinkage(X)
#' estimatedCov_QISshrink = cov_quadratic_inverse_shrinkage(X)
#' 
#' precision_NLshrink = solve(estimatedCov_NLshrink)
#' precision_QISshrink = solve(estimatedCov_QISshrink)
#' 
#' FrobeniusLoss2(precision_MoorePenrose_Cent, Sigma = Sigma)
#' FrobeniusLoss2(precision_MoorePenrose_NoCent, Sigma = Sigma)
#' FrobeniusLoss2(precision_NLshrink, Sigma = Sigma, type = "precision")
#' FrobeniusLoss2(precision_QISshrink, Sigma = Sigma, type = "precision")
#' 
#' # We now use the true value of the precision matrix as a target for shrinkage
#' precision_MoorePenrose_Cent_trueSigma = 
#'   Moore_Penrose_target(X, centeredCov = TRUE, Pi0 = solve(Sigma))
#' precision_MoorePenrose_NoCent_trueSigma = 
#'   Moore_Penrose_target(X, centeredCov = FALSE, Pi0 = solve(Sigma))                                                        
#'                                                         
#' FrobeniusLoss2(precision_MoorePenrose_Cent_trueSigma, Sigma = Sigma)
#' FrobeniusLoss2(precision_MoorePenrose_NoCent_trueSigma, Sigma = Sigma)
#' # this is indeed much closer than before
#' 
#' 
#' @export
Moore_Penrose_target <- function(X, centeredCov = TRUE, Pi0 = NULL, verbose = 0)
{
  if (is.null(Pi0)) {
    result = Moore_Penrose_target_general(X = X, centeredCov = centeredCov,
                                          Pi0 = Pi0, verbose = verbose)
  } else {
    result = Moore_Penrose_target_identity(X = X, centeredCov = centeredCov,
                                           verbose = verbose)
  }
  
  return (result)
}

Moore_Penrose_target_general <- function(X, Pi0 = NULL, centeredCov, verbose = 0)
{
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (is.null(Pi0)){
    Pi0 <- Ip
  } else if (nrow(Pi0) != p || ncol(Pi0) != p){
    stop("'Pi0' should be a 'p' by 'p' matrix.")
  }
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(X = X, centeredCov = centeredCov))
  
  
  ##### shrinkage MP
  trS1 <- tr(iS_MP) / p
  trS2 <- tr(iS_MP %*% iS_MP) / p
  trS3 <- tr(iS_MP %*% iS_MP %*% iS_MP) / p
  
  hv0 <- c_n * trS1
  ihv0 <- 1 / hv0
  ihv0_2 <- ihv0^2
  ihv0_3 <- ihv0^3
  
  h2 <- 1 / (trS2 * c_n)
  h3 <- trS3 / (trS2^3 * c_n^2)
  
  trS1Pi0 <- tr(iS_MP %*% Pi0) / p
  trS2Pi0 <- tr(iS_MP %*% iS_MP %*% Pi0) / p
  
  d0Pi0 <- tr( (Ip - S %*% iS_MP) %*% Pi0) / p
  # If Pi0 = Ip, then d0Pi0 should be equal to (c_n - 1) / c_n
  # This is the same as tr(Ip/p) - tr(S %*% iS_MP) / p
  # which is the same as 
  # tr(Ip/p) - tr(Ytilde %*% solve(t(Ytilde) %*% Ytilde) %*% t(Ytilde) ) / p
  # in the non-centered case.
  # 1 - tr(Ytilde %*% solve(t(Ytilde) %*% Ytilde) %*% t(Ytilde) ) / p
  # = 1 - tr(t(Ytilde) %*% Ytilde %*% solve(t(Ytilde) %*% Ytilde)  ) / p
  # = 1 - tr(diag(n)  ) / p
  #
  # But actually this is not correct, see:
  # tr(t(Ytilde) %*% Ytilde %*% solve(t(Ytilde) %*% Ytilde)  ) / p
  # tr(diag(n)  ) / p
  # This is because 
  # tr(t(Ytilde) %*% Ytilde %*% solve(t(Ytilde) %*% Ytilde)  ) / p = (n-1)/p
  #
  # So d0Pi0 
  # = 1 - (n-1) / p
  # = (p - (n-1)) / p
  # = (p/n - (1 - 1/n)) / (p/n)
  # = (c_n - (1 - 1/n)) / c_n
  #
  # while (c_n - 1) / c_n = (p / n - 1) / (p / n) = (p - n) / p
  # which is different.
  
  
  d1 <- trS1 / (trS2 * c_n)
  d1Pi0 <- trS1Pi0 / (trS2 * c_n)
  d2 <- (trS1 * trS3 - trS2 * trS2) / (trS2^3 * c_n^2)
  d2Pi0 <- (trS1Pi0 * trS3 - trS2 * trS2Pi0) / (trS2^3 * c_n^2)
  
  q1    <- tr(S) / p
  q1Pi0 <- tr(S %*% Pi0) / p
  q1Pi02 <- tr(S %*% Pi0 %*% Pi0) / p
  q2 <- tr(S %*% S) / p - c_n * q1 * q1
  q2Pi0 <- tr(S %*% S %*% Pi0) / p - c_n * q1 * q1Pi0
  q2Pi02 <- tr(S %*% S %*% Pi0 %*% Pi0) / p - c_n * q1 * q1Pi02
  
  d1Sig  <- ihv0 * (ihv0 / c_n - d1)
  # if (centeredCov){
  #   d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 * (1 - 1/n) /c_n)
  # } else {
  #   d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 / c_n)
  # }
  d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 / c_n)
  d1Sig2Pi0 <- ihv0_2 * (q1Pi0 + d1Pi0) - 2 * ihv0^3 * (tr(Pi0) / p - d0Pi0)
  
  d2Sig2 <- ihv0 * d1Sig2 - ihv0_2 * (d1Sig - d2)
  
  num_alpha_MP <- d1Sig * q2Pi0 - d1Sig2Pi0 * q1Pi0
  num_beta_MP <- -(d2Sig2 - d1Sig2 * h3 / h2) * q1Pi0 / h2 - d1Sig * d1Sig2Pi0 / h2
  den_MP   <- -(d2Sig2 - d1Sig2 * h3 / h2) * q2Pi02 / h2 - d1Sig2Pi0^2 / h2
  
  if (verbose > 0){
    cat("num_alpha_MP = ", num_alpha_MP, "\n")
    cat("num_beta_MP = ", num_beta_MP, "\n")
    cat("den_MP = ", den_MP, "\n")
  }
  
  alpha    <- num_alpha_MP / den_MP
  beta     <- num_beta_MP / den_MP
  iS_ShMP  <- alpha * iS_MP + beta * Pi0
  
  result = list(
    estimated_precision_matrix = iS_ShMP
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}


Moore_Penrose_target_identity <- function (X, centeredCov = TRUE, verbose = 0)
{
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(X = X, centeredCov = centeredCov))
  
  ##### shrinkage MP
  trS1<-sum(diag(iS_MP))/p
  trS2<-sum(diag(iS_MP%*%iS_MP))/p
  trS3<-sum(diag(iS_MP%*%iS_MP%*%iS_MP))/p
  
  hv0<-c_n*trS1
  ihv0<-1/hv0
  ihv0_2<-ihv0^2
  ihv0_3<-ihv0^3
  
  h2<-1/trS2/c_n
  h3<-trS3/(trS2^3)/c_n^2
  
  d1<-trS1/trS2/c_n
  d2<-(trS1*trS3-trS2^2)/(trS2^3)/c_n^2
  
  q1<-sum(diag(S))/p
  q2<-sum(diag(S%*%S))/p-c_n*(q1^2)
  
  d1Sig<-ihv0*(ihv0/c_n-d1)
  
  # if (centeredCov){
  #   d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 * (1 - 1/n) /c_n)
  # } else {
  #   d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 / c_n)
  # }
  d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 / c_n)
  d2Sig2<-ihv0*d1Sig2-ihv0_2*(d1Sig-d2)
  
  num_alpha_MP<-d1Sig*q2-d1Sig2*q1
  num_beta_MP<--(d2Sig2-d1Sig2*h3/h2)*q1/h2-d1Sig*d1Sig2/h2
  den_MP<--(d2Sig2-d1Sig2*h3/h2)*q2/h2-d1Sig2^2/h2
  
  if (verbose > 0){
    cat("num_alpha_MP = ", num_alpha_MP, "\n")
    cat("num_beta_MP = ", num_beta_MP, "\n")
    cat("den_MP = ", den_MP, "\n")
  }
  
  ha_MP <- num_alpha_MP / den_MP
  hb_MP <- num_beta_MP / den_MP
  iS_ShMP<-ha_MP * iS_MP + hb_MP * Ip
  
  result = list(
    estimated_precision_matrix = iS_ShMP
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return(result)
}



