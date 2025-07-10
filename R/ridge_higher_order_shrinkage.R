

#' Estimation of the derivative of the function v()
#' 
#' @param S_t_inverse_pow_jp1 inverse of the regularized covariance matrix
#' to the power of j+1
#' @param t regularization parameter
#' @param c_n p/n parameter (potentially adjusted for centering)
#' @param j order of the derivative of the function \eqn{v()}
#' 
#' @noRd
v_hat_j_of_t <- function(t, j, S_t_inverse_pow_jp1, c_n)
{
  p = nrow(S_t_inverse_pow_jp1)
  
  result = (-1)^j * factorial(j) * c_n * 
    ( tr(S_t_inverse_pow_jp1) / p - (c_n - 1)/c_n * t^(-(j+1) ) )
  
  # cat("j = ", j, "; trace = ", tr(S_t_inverse_pow_jp1) / p, "\n")
  
  return(result)
}


#' @param v_hat_until_j vector of size j with the first derivative of v until
#' the j-th derivative of v
#' 
#' @param h_hat_until_j vector of size j, with NA value at the first entry
#' and \eqn{\hat{h}_{k}} at the k-th entry
#' 
#' @param j value of j, has to be at least 2
#' 
#' @noRd
h_hat_jp1_t <- function(h_hat_until_j, v_hat_until_j, j)
{
  if (length(h_hat_until_j) != j){
    stop("length(h_hat_until_j) should be equal to j")
  }
  if (length(v_hat_until_j) != j){
    stop("length(v_hat_until_j) should be equal to j")
  }
  
  ugly_sum <- 0
  
  if (j > 2){
    # This is the sum from to 2 to j-1, so if j == 2, then the sum is empty
    # so we compute the sum only if j > 2.
    for (k in 2:(j-1))
    {
      ugly_sum <- ugly_sum +
        (-1)^(k) * factorial(k) * h_hat_until_j[k+1] *
        kStatistics::e_eBellPol(j, k, c(v_hat_until_j[1:(j - k + 1)], rep(0, k - 1)) )
    }
  }
  
  numerator = v_hat_until_j[j] + v_hat_until_j[1] * ugly_sum
  
  denominator = v_hat_until_j[1]^(j + 1) * (-1)^(j + 1) * factorial(j)
  
  result = numerator / denominator
  
  return (result)
}


#' @param h_hat_kmaxp1_t vector of size kmax, with the entries
#' k = 2 until k = kmax + 1
#' 
#' @noRd
compute_d_kl <- function(v_0_t, c_n, kmax, h_hat_kmaxp1_t, t, q1)
{
  if (length(h_hat_kmaxp1_t) != (kmax + 1 - 2 + 1)){
    stop("The length of h_hat_kmaxp1_t is ", length(h_hat_kmaxp1_t),
         "but it is supposed to be: ", (kmax + 1 - 2 + 1) ,
         "because it goes from 2 to ", kmax + 1, ".")
  }
  
  d_01 = 1 / (c_n * v_0_t) - t / c_n
  d_k1 = 1 / (c_n * v_0_t^(1 + (1:kmax)) ) - h_hat_kmaxp1_t / c_n
  
  d_02 = (q1 - d_01) / v_0_t
  
  d_k2 = rep(NA, kmax)
  for (k in 1:kmax){
    
    # Previous value of d_{k-1, 2} to be used in the recursion
    if (k == 1){
      d_km1_2 = d_02
    } else {
      d_km1_2 = d_k2[k - 1]
    }
    
    d_k2[k] <- (d_km1_2 - d_k1[k]) / v_0_t
  }
  
  result = cbind(c(d_01, d_k1), c(d_02, d_k2))
  colnames(result) <- c("l = 1", "l = 2")
  rownames(result) <- paste0("k = ", 0:(nrow(result) - 1) )
  
  return (result)
}



# Compute the matrix M for the higher-order shrinkage
# the output should be of size m+1
#
compute_M_t <- function(m, c_n, S_t_inverse, q1, q2, t, verbose)
{
  v_0_t <- v_hat_j_of_t(t = t, j = 0, S_t_inverse_pow_jp1 = S_t_inverse, c_n = c_n)
  
  v = rep(NA, 2 * m - 1)
  
  S_t_inverse_pow_jp1 = S_t_inverse
  
  for (j in 1:(2*m - 1)){
    S_t_inverse_pow_jp1 = S_t_inverse %*% S_t_inverse_pow_jp1
    v[j] <- v_hat_j_of_t(t = t, j = j,
                         S_t_inverse_pow_jp1 = S_t_inverse_pow_jp1,
                         c_n = c_n)
  }
  if (verbose){
    cat("Estimation of v(t0) = ", v_0_t, "\n\n")
    cat("Estimation of the derivatives of v:\n")
    print(v)
    cat("\n")
  }
  
  h <- rep(NA, 2 * m)
  h[2] = - 1 / v[1]
  
  if (2 <= 2 * m - 1){
    for (j in 2:(2 * m - 1)){
      h[j + 1] = h_hat_jp1_t(h_hat_until_j = h[1:j],
                             v_hat_until_j = v[1:j],
                             j = j)
    }
  }
  
  if (verbose){
    cat("Estimation of h:\n")
    print(h)
    cat("\n")
  }
  
  d <- compute_d_kl(v_0_t = v_0_t, c_n = c_n, kmax = 2 * m - 1,
                    h_hat_kmaxp1_t = h[2:(2 * m)],
                    t = t, q1 = q1)
  
  if (nrow(d) != 2 * m){
    stop("d should be of the right size.")
  }
  
  if (verbose){
    cat("Estimation of d[k,l]:\n")
    print(d)
    cat("\n")
  }
  
  s <- matrix(nrow = 2 * m, ncol = 2)
  for (index in 1:(2 * m)){
    j = index - 1
    # so that index = j + 1
    if (verbose){
      cat("We are now computing s[", index, ",], corresponding to j =", j, "\n",
          sep = "")
    }
    
    # Remember that d[1,l] corresponds in the paper to d_{0,l}
    # because indexing starts at 1.
    s[index, 1] = t^(-j - 1) * d[1, 1]
    s[index, 2] = t^(-j - 1) * d[1, 2]
    
    if (verbose){
      cat("d[0, 2] = ", d[1, 2], "\n")
      cat("s[", index, ", 2] =", s[index, 2], "\n", sep = "")
    }
    
    if (j >= 1){
      # Otherwise the sum is empty
      for (i in 1:j){
        for (k in 1:i){
          if (verbose){
            cat("This is term i =", i, ", k =", k, "\n")
          }
          Bell_polynomial = 
            kStatistics::e_eBellPol(i, k, c(v[1:(i - k + 1)], rep(0, k - 1)) )
          
          multiplicative_factor = t^(- (j - i) - 1) * (-1)^(i + k) * 
            factorial(k) / factorial(i) * Bell_polynomial
          
          if (verbose){
            cat("Bell_polynomial = ", Bell_polynomial, "\n")
          }
          
          for (l in 1:2){
            # Remember that d[k + 1, l] corresponds in the paper to d_{k,l}
            additional_term = multiplicative_factor * d[k + 1, l]
            
            if (verbose && (l == 2)){
              cat("l = ", l, ", d_{k,l} = ", d[k + 1, l], "\n")
              cat("l = ", l, ", additional_term = ", additional_term, "\n")
            }
            
            s[index, l] = s[index, l] + additional_term
          }
        }
      }
    }
  }
  
  if (verbose){
    cat("Estimation of s:\n")
    print(s)
    cat("\n")
  }
  
  # Computation of M  ==========================================================
  s2 = c(q2, s[, 2])
  
  M <- s2[1:(m+1)]
  for (j in 2:(m+1))
  {
    M <- cbind(M, s2[j:(m+j)])
  }
  
  
  # Computation of hm  =========================================================
  hm = c(q1, s[1:m, 1])
  
  if (verbose){
    cat("Estimation of hm:\n")
    print(hm)
    cat("\n")
  }
  
  return(list(M = M, hm = hm, v = v))
}



#' Ridge higher order shrinkage
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
#' precision_MoorePenrose_Cent = Moore_Penrose_shrinkage(Y = t(X), centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = Moore_Penrose_shrinkage(Y = t(X), centeredCov = FALSE)
#'
#' FrobeniusLoss2(precision_MoorePenrose_Cent, Sigma = Sigma)
#' FrobeniusLoss2(precision_MoorePenrose_NoCent, Sigma = Sigma)
#' 
#' for (m in 1:5){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       ridge_higher_order_shrinkage(Y = t(X), m = m, centeredCov = TRUE, t = 10)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       ridge_higher_order_shrinkage(Y = t(X), m = m, centeredCov = FALSE, t = 10)
#'       
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_Cent, Sigma = Sigma))
#'   
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_NoCent, Sigma = Sigma))
#' }
#' 
#' 
#' precision_higher_order_shrinkage_Cent = 
#'   ridge_higher_order_shrinkage(Y = t(X), m = 1, centeredCov = TRUE, t = 100)
#' 
#' precision_target_identity_semioptimal_Cent = 
#'   ridge_target_identity_semioptimal(Y = t(X), centeredCov = TRUE, t = 100)
#'       
#' precision_target_general_semioptimal_Cent = 
#'   ridge_target_general_semioptimal(Y = t(X), centeredCov = TRUE, t = 100, Pi0 = diag(p))
#'
#' precision_higher_order_shrinkage_Cent$alpha
#' precision_target_identity_semioptimal_Cent$beta_optimal
#' precision_target_identity_semioptimal_Cent$alpha_optimal
#' 
#' precision_target_general_semioptimal_Cent$beta_optimal
#' precision_target_general_semioptimal_Cent$alpha_optimal
#' 
#' precision_higher_order_shrinkage_Cent$M
#' precision_target_identity_semioptimal_Cent$M
#' 
#' precision_higher_order_shrinkage_Cent$hm
#' precision_target_identity_semioptimal_Cent$hm
#' 
#' FrobeniusLoss2(precision_higher_order_shrinkage_Cent, Sigma = Sigma)
#' FrobeniusLoss2(precision_target_identity_semioptimal_Cent, Sigma = Sigma)
#' 
#' 
#' @export
#' 
ridge_higher_order_shrinkage <- function(Y, m, centeredCov, t, verbose = 2)
{
  if (verbose){
    cat("Starting `ridge_higher_order_shrinkage`...\n")
  }
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  
  if (verbose){
    cat("*  n = ", n, "\n")
    cat("*  p = ", p, "\n")
    cat("*  t = ", t, "\n")
    cat("*  m = ", m, "\n")
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
    
    c_n = p / (n-1)
    
  } else {
    if (verbose){
      cat("*  non-centered case\n")
    }
    
    S <- Y %*% t(Y)/n
    
    c_n = p / n
  }
  if (verbose){
    cat("*  c_n = ", c_n, "\n\n")
  }
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  S_t <- S + t * diag(nrow = p)
  
  S_t_inverse <- solve(S_t)
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  if (verbose){
    cat("Starting values: \n")
    cat("*  q1 = ", q1, "\n")
    cat("*  q2 = ", q2, "\n\n")
  }
  
  estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                           S_t_inverse = S_t_inverse,
                           t = t, verbose = verbose)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  if (verbose){
    cat("Optimal alpha: \n")
    print(alpha)
  }
  
  estimated_precision_matrix = alpha[1] * Ip
  power_S_t_inverse = Ip
  
  for (k in 1:m){
    power_S_t_inverse = power_S_t_inverse %*% S_t_inverse
    estimated_precision_matrix = estimated_precision_matrix + 
      alpha[k + 1] * power_S_t_inverse
  }
  
  result = list(
    estimated_precision_matrix = estimated_precision_matrix,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

