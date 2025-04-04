
#' First-order shrinkage of the Moore-Penrose inverse towards a general target
#'
#'
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param Pi0 prior of the precision matrix. This a `p` by `p` matrix, used as
#' a target for the shrinkage. Default value is the identity matrix of size `p`.
#' As an advice, it should be a symmetric positive-definite matrix, but this is
#' not checked for.
#' 
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix, i.e. the inverse of the covariance matrix).
#' It is (asymptotically) optimal for the loss
#' \eqn{Loss(EstimatorPi) := \| EstimatorPi * \Sigma - I \|_F^2}.
#' 
#' 
#' @examples
#' p = 500
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
#' precision_MoorePenrose_Cent = Moore_Penrose_shrinkage(t(X), centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = Moore_Penrose_shrinkage(t(X), centeredCov = FALSE)
#' precision_MoorePenrose_toIPCent = Moore_Penrose_shrinkage_toIP(t(X), centeredCov = TRUE)
#' precision_MoorePenrose_toIPNoCent = Moore_Penrose_shrinkage_toIP(t(X), centeredCov = FALSE)
#' 
#' precisionTrue = solve(Sigma)
#' 
#' estimatedCov_NLshrink = analytical_NL_shrinkage(t(X))
#' estimatedCov_QISshrink = quadratic_inverse_shrinkage(X)
#' 
#' precision_NLshrink = solve(estimatedCov_NLshrink)
#' precision_QISshrink = solve(estimatedCov_QISshrink)
#' 
#' FrobeniusNorm2 <- function(M){sum(diag(M %*% t(M)))}
#' FrobeniusNorm2(precision_MoorePenrose_Cent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_MoorePenrose_NoCent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_MoorePenrose_toIPCent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_MoorePenrose_toIPNoCent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_NLshrink     %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_QISshrink    %*% Sigma - diag(p) ) / p
#' 
#' # We now use the true value of the precision matrix as a target for shrinkage
#' precision_MoorePenrose_Cent = Moore_Penrose_shrinkage(t(X), centeredCov = TRUE, Pi0 = solve(Sigma))
#' precision_MoorePenrose_NoCent = Moore_Penrose_shrinkage(t(X), centeredCov = FALSE, Pi0 = solve(Sigma))
#' FrobeniusNorm2(precision_MoorePenrose_Cent %*% Sigma - diag(p) ) / p
#' FrobeniusNorm2(precision_MoorePenrose_NoCent %*% Sigma - diag(p) ) / p
#' # this is indeed much closer than before
#' 
#' 
Moore_Penrose_shrinkage <- function(Y, Pi0 = NULL, centeredCov)
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = p / n
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (is.null(Pi0)){
    Pi0 <- Ip
  } else if (nrow(Pi0) != p || ncol(Pi0) != p){
    stop("'Pi0' should be a 'p' by 'p' matrix.")
  }
  
  # Inverse companion covariance
  iY <- solve(t(Y) %*% Y / n)
  
  if (centeredCov){
    S <- cov(t(Y))
  } else {
    S <- Y %*% t(Y)/n
  }
  eigenS <- eigen(S)
  U <- eigenS$vectors
  
  ##### Moore-Penrose inverse
  iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
  D_MP <- diag(eigen(iS_MP)$values)
  
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
  d1Sig2 <- ihv0_2 * (q1 + d1 - 2 * ihv0 / c_n)
  d1Sig2Pi0 <- ihv0_2 * (q1Pi0 + d1Pi0) - 2 * ihv0^3 * (tr(Pi0) / p - d0Pi0)
  
  d2Sig2 <- ihv0 * d1Sig2 - ihv0_2 * (d1Sig - d2)
  
  num_alpha_MP <- d1Sig * q2 - d1Sig2Pi0 * q1Pi0
  num_beta_MP <- -(d2Sig2 - d1Sig2 * h3 / h2) * q1Pi0 / h2 - d1Sig * d1Sig2Pi0 / h2
  den_MP   <- -(d2Sig2 - d1Sig2 * h3 / h2) * q2Pi02 / h2 - d1Sig2Pi0^2 / h2
  
  alpha    <- num_alpha_MP / den_MP
  beta     <- num_beta_MP / den_MP
  iS_ShMP  <- alpha * iS_MP + beta * Pi0
  
  return (iS_ShMP)
}



Moore_Penrose_shrinkage_toIP <- function (Y, centeredCov)
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = p / n
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Inverse companion covariance
  iY <- solve(t(Y) %*% Y / n)
  
  if (centeredCov){
    S <- cov(t(Y))
  } else {
    S <- Y %*% t(Y)/n
  }
  
  eigenS <- eigen(S)
  U <- eigenS$vectors
  
  ##### Moore-Penrose inverse
  iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
  D_MP <- diag(eigen(iS_MP)$values)
  
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
  d1Sig2<-ihv0_2*(q1+d1-2*ihv0/c_n)
  d2Sig2<-ihv0*d1Sig2-ihv0_2*(d1Sig-d2)
  
  num_a_MP<-d1Sig*q2-d1Sig2*q1
  num_b_MP<--(d2Sig2-d1Sig2*h3/h2)*q1/h2-d1Sig*d1Sig2/h2
  den_MP<--(d2Sig2-d1Sig2*h3/h2)*q2/h2-d1Sig2^2/h2
  ha_MP<-num_a_MP/den_MP
  hb_MP<-num_b_MP/den_MP
  iS_ShMP<-ha_MP*iS_MP+hb_MP*Ip
  
  return(iS_ShMP)
}



