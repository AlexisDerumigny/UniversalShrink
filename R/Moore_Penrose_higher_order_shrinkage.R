

#' Estimate the derivatives of the function v at t = o
#' 
#' @param m parameter describing the derivatives to take. We want to estimate
#' all the derivatives from i = 1 to i = m
#' 
#' @param c_n ratio p/n, eventually adjusted for centering
#' 
#' @param p number of variables
#' 
#' @param Splus Moore-Penrose inverse of the sample covariance matrix
#' 
#' @returns a vector of size \code{m} containing
#' \eqn{\widehat{v'}(0), \dots, \widehat{v^{(m)}}(0)}.
#' 
#' @noRd
estimator_vhat_derivative_t0 <- function(m, c_n, p, Splus){
  v <- rep(NA, m)
  
  for (i in 1:m)
  {
    v[i] <- (-1)^i * factorial(i) * c_n * (1/p) * tr(Splus^(i+1))
  }
  
  return (v)
} 

#' This function estimates \eqn{\widehat{v(0)}}
#'
#' @noRd
estimator_vhat_t0 <- function(c_n, p, Splus){
  v <- c_n * (1/p) * tr(Splus)
  
  return (v)
}

# Compute the matrix M for the higher-order shrinkage
compute_M <- function(m, n, p, ihv0, D_MP, q1, q2, h2, h3, hv0, centeredCov)
{
  if (centeredCov){
    c_n = p / (n-1)
  } else {
    c_n = p/n
  }
  
  h <- rep(NA, 2*m+1)
  d <- matrix(NA, nrow = 2 * m, ncol = 2)
  
  h[2] <- h2
  h[3] <- h3
  
  v <- estimator_vhat_derivative_t0(m = 2 * m, c_n, p = p, Splus = D_MP)
  
  if (m > 1)
  { 
    for (i in 3:(2*m))
    {
      sh<- 0
      for (k in 2:(i-1))
      {
        sh <- sh + (-1)^k * factorial(k) * h[k+1] * 
          kStatistics::e_eBellPol(i, k, c(v[1:(i-k+1)], rep(0,k-1)))
      }
      h[i+1] <- (v[i] + v[1] * sh) / ((v[1])^(i+1) * (-1)^(i+1) * factorial(i))
    }
  }
  
  for (k in 1:(2*m))
  {
    for (l in 1:2)
    {
      if (l == 1){
        d[k,l] <- 1 / c_n * (1 / (hv0^(k+1)) - h[k+1])
      } else if (l == 2 && k == 1) {
        # if (centeredCov){
        #   d[k,l] <- ihv0 * (ihv0 * (q1 - ihv0 * ( (1 - 1/n) / c_n ) ) - d[1,1] )
        # } else {
        #   d[k,l] <- ihv0 * (ihv0 * (q1 - ihv0 * ( 1 / c_n ) ) - d[1,1] )
        # }
        
        d[k,l] <- ihv0 * (ihv0 * (q1 - ihv0 * ( 1 / c_n ) ) - d[1,1] )
      } else if (l == 2 && k > 1) {
        d[k,l] <- ihv0 * (d[k-1,2] - d[k,1])
      }
    }
  }
  
  s2 <- rep(NA , 2 * m + 1)
  s2[1] <- q2
  
  for (j in 1:(2*m))
  {
    s2[j+1] <- 0
    for (k in 1:j)
    {
      s2[j+1] <- s2[j+1] +
        (-1)^(j + k + 1) * factorial(k) / factorial(j) * d[k,2] * 
        kStatistics::e_eBellPol(j, k, c(v[1:(j - k + 1)], rep(0, k - 1)))
    }
  }
  
  M <- s2[1:(m+1)]
  for (j in 2:(m+1))
  {
    M <- cbind(M, s2[j:(m+j)])
  }
  
  M <- apply(M, c(1,2), Re)
  
  
  # Computation of hm
  
  hm <- rep(NA, m + 1)
  
  hm[1] <- q1
  
  for (j in 1:m)
  {
    hm[j+1]<-0
    for (k in 1:j)
    {
      hm[j+1] <- hm[j+1] +
        (-1)^(j + k + 1) * factorial(k) / factorial(j) * d[k,1] *
        kStatistics::e_eBellPol(j, k, c(v[1:(j - k + 1)], rep(0, k - 1)))
    }
  }
  hm <- Re(hm)
  
  return(list(M = M, hm = hm, v = v))
}



#' Moore-Penrose higher order shrinkage of the precision matrix
#' 
#' This function computes an estimator of the precision matrix by using the
#' polynomial
#' \deqn{
#' \mathbf{S}_{n;HOS}^+=\hat{\alpha}_0^+\mathbf{I}_p 
#' +\sum_{j=1}^m\hat{\alpha}^+_j(\mathbf{S}^+_n)^j,
#' }
#' where \eqn{\hat{\boldsymbol{\alpha}}^{+}(m)=(\hat{\alpha}_0^+,\hat{\alpha}_1^+,
#' \ldots,\hat{\alpha}_m^+)^\top} given by
#'\deqn{
#'\hat{\boldsymbol{\alpha}}^{+}(m)=\widehat{\mathbf{M}}^+(m)^{-1} \hat{\mathbf{m}}^+(m)
#'}  with
#' \deqn{
#'\hat{\mathbf{m}}^+(m)=
#'  \begin{pmatrix}
#'\hat{q}_1\\
#'\hat{s}_{1,1}\\
#'\vdots\\
#'\hat{s}_{m,1}
#'\end{pmatrix}
#'\quad \text{and}\quad
#'\widehat{\mathbf{M}}^+(m)=\begin{pmatrix}
#'\hat{q}_2 & \hat{s}_{1,2}  & \ldots & \hat{s}_{m,2} \\
#'\hat{s}_{1,2}& \hat{s}_{2,2}  & \ldots & \hat{s}_{m+1,2} \\
#'\vdots&\vdots&\ddots&\vdots\\
#'\hat{s}_{m,2}& \hat{s}_{m+1,2} & \ldots & \hat{s}_{2m,2}
#'\end{pmatrix},
#'} where
#' \eqn{\mathbf{S}^+_n} is the Moore-Penrose inverse of the sample
#' covariance matrix and \eqn{\mathbf{I}_p} is the identity matrix of size \eqn{p}.
#' The details on the computation of the terms \eqn{\hat q_1}, \eqn{\hat q_2} and
#'  \eqn{s_{i,j}} are given in 
#' Theorem 2.5 of Bodnar and Parolya (2025).
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix, i.e. the inverse of the covariance matrix).
#' 
#' @references
#' Nestor Parolya & Taras Bodnar (2025).
#' Higher-order nonlinear shrinkage estimator of large-dimensional precision matrix.
#' \doi{10.1090/tpms/1239}
#' 
#' @examples
#' 
#' n = 10
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
#' precision_MoorePenrose_Cent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = FALSE)
#'
#' print(FrobeniusLoss2(precision_MoorePenrose_Cent, Sigma))
#' print(FrobeniusLoss2(precision_MoorePenrose_NoCent, Sigma))
#' 
#' for (m in 1:5){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       Moore_Penrose_higher_order_shrinkage(Y = t(X), m = m, centeredCov = TRUE)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       Moore_Penrose_higher_order_shrinkage(Y = t(X), m = m, centeredCov = FALSE)
#'       
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_Cent, Sigma))
#'   
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_NoCent, Sigma))
#' }
#' 
#' 
#' 
#' @export
#' 
Moore_Penrose_higher_order_shrinkage <- function(Y, m, centeredCov)
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
  iS_MP = as.matrix(Moore_Penrose(Y, centeredCov = centeredCov))
  
  D_MP <- diag(eigen(iS_MP)$values)
  
  trS1 <- tr(iS_MP) / p
  trS2 <- tr(iS_MP %*% iS_MP) / p
  trS3 <- tr(iS_MP %*% iS_MP %*% iS_MP) / p
  
  hv0 <- estimator_vhat_t0(c_n = c_n, p = p, Splus = D_MP)
  ihv0 <- 1 / hv0
  ihv0_2 <- ihv0^2
  ihv0_3 <- ihv0^3
  
  h2 <- 1 / (trS2 * c_n)
  h3 <- trS3 / (trS2^3 * c_n^2)
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  estimatedM = compute_M(m = m, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
                         q1 = q1, q2 = q2, h2 = h2, h3 = h3, hv0 = hv0,
                         centeredCov = centeredCov)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  
  result = alpha[1] * Ip
  power_isMP = Ip
  
  for (k in 1:m){
    power_isMP = power_isMP %*% iS_MP
    result = result + alpha[k + 1] * power_isMP
  }
  
  result = list(
    estimated_precision_matrix = result,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

