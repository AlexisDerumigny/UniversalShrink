

#' Estimate the derivatives of the function v at t = o
#' 
#' @param m parameter describing the derivatives to take. We want to estimate
#' all the derivatives from i = 1 to i = m
#' 
#' @param c_n ratio p/n, eventually adjusted for centering
#' 
#' @param p number of variables
#' 
#' @param D_MP diagonal matrix containing on the diagonal the eigenvalues of the
#' Moore-Penrose inverse of the sample covariance matrix.
#' 
#' @returns a vector of size \code{m} containing
#' \eqn{\widehat{v'}(0), \dots, \widehat{v^{(m)}}(0)}.
#' 
#' @noRd
estimator_vhat_derivative_t0 <- function(m, c_n, p, D_MP){
  v <- rep(NA, m)
  
  for (i in 1:m)
  {
    v[i] <- (-1)^i * factorial(i) * c_n * (1/p) * sum(diag(D_MP)^(i+1))
  }
  
  return (v)
}

estimator_w_hat_derivative_t0 <- function(m, c_n, p, D_MP){
  w <- rep(NA, m)
  
  for (i in 1:m)
  {
    w[i] <- (-1)^(i-1) * factorial(i) * c_n * (1/p) * sum(diag(D_MP)^(i))
  }
  return (w)
} 

#' This function estimates \eqn{\widehat{v(0)}}
#'
#' @noRd
estimator_vhat_t0 <- function(c_n, p, D_MP){
  v <- c_n * (1/p) * tr(D_MP)
  
  return (v)
}

estimator_w_hat_tilde_derivative_t0 <- function(m, c_n, w_hat){
  wtilde <- rep(NA, m)
  wtilde[1] <- 1 / (1 - c_n)
  
  for (j in 2:m)
  {
    wtilde[j] <- 0
    for (k in 1:(j-1)) {
      Bell_polynomial = 
        kStatistics::e_eBellPol(j - 1, k, c(w_hat[1:(j - k)], rep(0, k - 1)))
      
      additional_term = 
        ((-1)^k * factorial(k) / (1 - c_n)^(k + 1)) * Bell_polynomial
      
      wtilde[j] <- wtilde[j] + additional_term
    }
    wtilde[j] <- j * wtilde[j]
  }
  return (wtilde)
} 


estimator_d_hat_tilde_p_small <- function(m, c_n, w_hat_tilde, Bell_polynomials){
  d = rep(NA, m)
  for (j in 1:m){
    second_term = 0
    if (j > 1){
      for (k in 1:(j - 1)) {
        # Bell_polynomial = 
        #   kStatistics::e_eBellPol(j, k, c(w_hat_tilde[1:(j - k + 1)], rep(0, k - 1)))
        Bell_polynomial = Bell_polynomials[j, k]
        
        new_term = (-1)^k * factorial(k) * d[k] * Bell_polynomial
        second_term = second_term + new_term
      }
      second_term = second_term * c_n
    }
    num = w_hat_tilde[j] + second_term
    den = c_n * (-1)^(j + 1) * factorial(j) * (1 - c_n)^(- j)
    d[j] <- num / den
  }
  return (d)
}


# Compute the matrix M for the higher-order shrinkage in the p > n regime.
compute_M_MoorePenrose_plarge <- function(
    m, n, p, c_n, ihv0, q1, q2, h2, h3, hv0, centeredCov, D_MP, method_invM, verbose)
{
  if (verbose > 0){
    cat("Starting compute_M_MoorePenrose_plarge...\n")
  }
  
  v <- estimator_vhat_derivative_t0(m = 2 * m, c_n, p = p, D_MP = D_MP)
  
  # Bell_polynomials = matrix(nrow = 2 * m, ncol = 2 * m)
  # for (j in 1:(2*m))
  # {
  #   for (k in 1:j)
  #   {
  #     Bell_polynomials[j, k] = 
  #       kStatistics::e_eBellPol(j, k, c(v[1:(j - k + 1)], rep(0, k - 1)))
  #   }
  # }
  Bell_polynomials = bellPolynomials(v, verbose = 0)
  # Removing the lines corresponding to n = 0 and k = 0
  Bell_polynomials = Bell_polynomials[-1, -1]
  
  h <- rep(NA, 2*m+1)
  h[2] <- h2
  h[3] <- h3
  
  if (m > 1)
  { 
    for (i in 3:(2*m))
    {
      sh<- 0
      for (k in 2:(i-1))
      {
        sh <- sh + (-1)^k * factorial(k) * h[k+1] * Bell_polynomials[i, k]
      }
      h[i+1] <- (v[i] + v[1] * sh) / ((v[1])^(i+1) * (-1)^(i+1) * factorial(i))
    }
  }
  
  d <- matrix(NA, nrow = 2 * m, ncol = 2)
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
      s2[j+1] <- s2[j+1] + (-1)^(j + k + 1) * factorial(k) / factorial(j) * 
        d[k,2] * Bell_polynomials[j, k]
    }
  }
  
  if (method_invM == "solve"){
    if (verbose > 0){
      cat("Numerical inversion of the matrix M...\n")
    }
    # M <- s2[1:(m+1)]
    # for (j in 2:(m+1))
    # {
    #   M <- cbind(M, s2[j:(m+j)])
    # }
    
    ## Fast generation of Hankel matrices
    id <- 1:(m+1) + rep(0:m, each = m + 1)
    
    M <- matrix(s2[id], nrow = m + 1)
    
    invM = solve(M)
    
    if (verbose > 1){
      cat("M = \n")
      print(M)
      cat("M^{-1} = \n")
      print(invM)
    }
    
  } else if (method_invM == "recursive"){
    if (verbose > 0){
      cat("Using the recursive formula to compute the inverse of the matrix M...\n")
    }
    # M is not computed but we still need to declare this variable because it
    # is used in the returned list.
    M = NULL
    
    # We avoid computing M and inverting it numerically. Here we compute the
    # inverse of the matrix M by using the recursive formula.
    invM = compute_M_inverse(m = m, all_tr0 = 1 / s2[1],
                             all_tr = s2[-1], verbose = verbose - 2)
    if (verbose > 1){
      cat("M^{-1} = \n")
      print(invM)
    }
  } else {
    stop("method_invM '", method_invM, "' unavailable. Possible choices are: ",
         "'solve' and 'recursive'.")
  }
  
  
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
        Bell_polynomials[j, k]
    }
  }
  hm <- Re(hm)
  if (verbose > 0){
    cat("hm = \n")
    print(hm)
  }
  
  alpha = invM %*% hm
  
  return(list(M = M, invM = invM, hm = hm, v = v, alpha = alpha))
}



# Compute the matrix M for the higher-order shrinkage in the p < n regime.
compute_M_MoorePenrose_psmall <- function(
    m, n, p, c_n, q1, centeredCov, D_MP, method_invM, verbose)
{
  if (verbose > 0){
    cat("Starting compute_M_MoorePenrose_psmall...\n")
  }
  w_hat <- estimator_w_hat_derivative_t0(
    m = 2 * m, c_n = c_n, p = p, D_MP = D_MP)
  
  w_hat_tilde <- estimator_w_hat_tilde_derivative_t0(
    m = 2 * m, c_n = c_n, w_hat = w_hat)
  
  
  # Bell_polynomials = matrix(nrow = 2 * m, ncol = 2 * m)
  # for (j in 1:(2*m))
  # {
  #   for (k in 1:j)
  #   {
  #     Bell_polynomials[j, k] = kStatistics::e_eBellPol(
  #     j, k, c(w_hat_tilde[1:(j - k + 1)], rep(0, k - 1) ) )
  #   }
  # }
  Bell_polynomials = bellPolynomials(w_hat_tilde, verbose = 0)
  # Removing the lines corresponding to n = 0 and k = 0
  Bell_polynomials = Bell_polynomials[-1, -1]
  
  
  d_hat_tilde_positive <- estimator_d_hat_tilde_p_small(
    m = 2 * m, c_n = c_n, w_hat_tilde = w_hat_tilde,
    Bell_polynomials = Bell_polynomials)
  
  d_hat_tilde_0 = 1
  d_hat_tilde_minus1 = q1
  
  # Function that represents an array from -1 to 2 * m
  d_hat_tilde = function(k){
    if (k == -1){
      return (d_hat_tilde_minus1)
    } else if (k == 0){
      return (d_hat_tilde_0)
    } else {
      return (d_hat_tilde_positive[k])
    }
  }
  
  s2 <- rep(NA , 2 * m + 1)
  
  s2[1] <- q1 / (1 - c_n)
  
  for (j in 1:(2*m))
  {
    s2[j+1] <- 0
    for (k in 1:j)
    {
      s2[j+1] <- s2[j+1] +
        (-1)^(j + k) * factorial(k) / factorial(j) * d_hat_tilde(k - 2) * 
        Bell_polynomials[j, k]
    }
  }
  
  if (method_invM == "solve"){
    if (verbose > 0){
      cat("Numerical inversion of the matrix M...\n")
    }
    # M <- s2[1:(m+1)]
    # for (j in 2:(m+1))
    # {
    #   M <- cbind(M, s2[j:(m+j)])
    # }
    
    ## Fast generation of Hankel matrices
    id <- 1:(m+1) + rep(0:m, each = m + 1)
    
    M <- matrix(s2[id], nrow = m + 1)
    
    invM = solve(M)
    
    if (verbose > 1){
      cat("M = \n")
      print(M)
      cat("M^{-1} = \n")
      print(invM)
    }
  } else if (method_invM == "recursive"){
    if (verbose > 0){
      cat("Using the recursive formula to compute the inverse of the matrix M...\n")
    }
    # M is not computed but we still need to declare this variable because it
    # is used in the returned list.
    M = NULL
    
    # We avoid computing M and inverting it numerically. Here we compute the
    # inverse of the matrix M by using the recursive formula.
    invM = compute_M_inverse(m = m, all_tr0 = 1 / s2[1],
                             all_tr = s2[-1], verbose = verbose - 2)
    if (verbose > 1){
      cat("M^{-1} = \n")
      print(invM)
    }
  } else {
    stop("method_invM '", method_invM, "' unavailable. Possible choices are: ",
         "'solve' and 'recursive'.")
  }
  
  
  # Computation of hm. This is s_{j,1}.
  
  hm <- rep(NA, m + 1)
  
  hm[1] <- q1
  hm[2] <- 1 / (1 - c_n)
  
  if (m >= 2){ # Indeed if m = 1 there is nothing to be done.
    for (j in 2:m)
    {
      hm[j+1]<-0
      for (k in 1:j)
      {
        hm[j+1] <- hm[j+1] +
          (-1)^(j + k) * factorial(k) / factorial(j) * d_hat_tilde(k - 1) *
          Bell_polynomials[j, k]
      }
    }
  }
  
  hm <- Re(hm)
  if (verbose > 0){
    cat("hm = \n")
    print(hm)
  }
  
  alpha = invM %*% hm
  
  return(list(M = M, invM = invM, hm = hm, alpha = alpha))
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
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' @param method_invM method for computing the inverse of the matrix M.
#' It can be \code{"solve"} (computing M and then inverting it) or
#' \code{"recursive"} (using a recursive formula, which is more stable
#' numerically for large \eqn{m}, and is therefore the default).
#' 
#' @inheritParams cov_with_centering
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
#' precision_MoorePenrose_Cent = Moore_Penrose_shrinkage(X, centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = Moore_Penrose_shrinkage(X, centeredCov = FALSE)
#'
#' print(LossFrobenius2(precision_MoorePenrose_Cent, Sigma))
#' print(LossFrobenius2(precision_MoorePenrose_NoCent, Sigma))
#' 
#' for (m in 1:3){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = TRUE)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = FALSE)
#'       
#'   print(LossFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   
#'   print(LossFrobenius2(precision_higher_order_shrinkage_NoCent, Sigma))
#' }
#' 
#' 
#' # Example for the case p < n
#' n = 10
#' p = n / 2
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
#' precision_MoorePenrose_Cent = Moore_Penrose_shrinkage(X, centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = Moore_Penrose_shrinkage(X, centeredCov = FALSE)
#'
#' print(LossFrobenius2(precision_MoorePenrose_Cent, Sigma))
#' print(LossFrobenius2(precision_MoorePenrose_NoCent, Sigma))
#' 
#' for (m in 1:3){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = TRUE)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       Moore_Penrose_higher_order_shrinkage(X, m = m, centeredCov = FALSE)
#'       
#'   print(LossFrobenius2(precision_higher_order_shrinkage_Cent, Sigma))
#'   
#'   print(LossFrobenius2(precision_higher_order_shrinkage_NoCent, Sigma))
#' }
#' 
#' @export
#' 
Moore_Penrose_higher_order_shrinkage <- function(
    X, m, centeredCov = TRUE, method_invM = "recursive", verbose = 0)
{
  call_ = match.call()
  # Get sizes of X
  n = nrow(X)
  p = ncol(X)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = X, centeredCov = centeredCov)
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- Moore_Penrose(X = X, centeredCov = centeredCov)
  iS_MP_ <- as.matrix(iS_MP)
  
  D_MP <- diag(eigen(iS_MP)$values)
  q1 <- tr(S) / p
  
  if (p < iS_MP$n_adjusted){
    estimatedM = compute_M_MoorePenrose_psmall(
      m = m, n = n, p = p, c_n = c_n, q1 = q1,
      centeredCov = centeredCov, D_MP = D_MP, method_invM = method_invM,
      verbose = verbose)
    
  } else if (p > iS_MP$n_adjusted){
    
    # trS1 <- tr(iS_MP_) / p
    trS2 <- tr(iS_MP_ %*% iS_MP_) / p
    trS3 <- tr(iS_MP_ %*% iS_MP_ %*% iS_MP_) / p
    
    hv0 <- estimator_vhat_t0(c_n = c_n, p = p, D_MP = D_MP)
    ihv0 <- 1 / hv0
    ihv0_2 <- ihv0^2
    ihv0_3 <- ihv0^3
    
    h2 <- 1 / (trS2 * c_n)
    h3 <- trS3 / (trS2^3 * c_n^2)
    q2 <- tr(S %*% S) / p - c_n * q1^2
    
    
    estimatedM = compute_M_MoorePenrose_plarge(
      m = m, n = n, p = p, c_n = c_n, ihv0 = ihv0,
      q1 = q1, q2 = q2, h2 = h2, h3 = h3, hv0 = hv0,
      centeredCov = centeredCov, D_MP = D_MP, method_invM = method_invM,
      verbose = verbose)
  } else {
    stop("This estimator is not defined for p = n - 1 in the centered case,",
         "and for p = n in the non-centered case. Here p = ", p,
         "and n = ", n, ".")
  }
  
  
  alpha = estimatedM$alpha
  
  result = alpha[1] * Ip
  power_isMP = Ip
  
  for (k in 1:m){
    power_isMP = power_isMP %*% iS_MP_
    result = result + alpha[k + 1] * power_isMP
  }
  
  result = list(
    estimated_precision_matrix = result,
    M = estimatedM$M,
    invM = estimatedM$invM,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v,
    n = n,
    p = p,
    centeredCov = centeredCov,
    method_invM = method_invM,
    method = "Moore-Penrose higher-order shrinkage",
    call = call_
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

