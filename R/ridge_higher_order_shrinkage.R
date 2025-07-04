

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
    
    d_k2[k] <- (d_km1_2 - d_k1[k])
  }
  
  result = cbind(c(d_01, d_k1), c(d_02, d_k2))
  
  return (result)
}



# Compute the matrix M for the higher-order shrinkage
# the output should be of size m+1
#
compute_M_t <- function(m, c_n, S_t_inverse, q1, q2, t)
{
  v_0_t <- v_hat_j_of_t(t = t, j = 0, S_t_inverse_pow_jp1 = S_t_inverse, c_n = c_n)
  
  v = rep(NA, 2 * m)
  
  S_t_inverse_pow_jp1 = S_t_inverse
  
  for (j in 1:(2*m)){
    S_t_inverse_pow_jp1 = S_t_inverse %*% S_t_inverse_pow_jp1
    v[j] <- v_hat_j_of_t(t = t, j = j,
                         S_t_inverse_pow_jp1 = S_t_inverse_pow_jp1,
                         c_n = c_n)
  }
  
  h <- rep(NA, 2*m+1)
  h[2] = - 1 / v[1]
  
  for (j in 2:(2*m)){
    h[j + 1] = h_hat_jp1_t(h_hat_until_j = h[1:j],
                           v_hat_until_j = v[1:j],
                           j = j)
  }
  
  
  d <- compute_d_kl(v_0_t = v_0_t, c_n = c_n, kmax = 2 * m,
                    h_hat_kmaxp1_t = h[2:(2*m+1)],
                    t = t, q1 = q1)
  
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
#' FrobeniusNorm2 <- function(M){sum(diag(M %*% t(M)))}
#' 
#' precision_MoorePenrose_Cent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = TRUE)
#' precision_MoorePenrose_NoCent = 
#'   Moore_Penrose_shrinkage(Y = t(X), centeredCov = FALSE)
#'
#' FrobeniusLoss2 <- function(M){FrobeniusNorm2(M %*% Sigma - diag(p) ) / p}
#'
#' print(FrobeniusLoss2(precision_MoorePenrose_Cent))
#' print(FrobeniusLoss2(precision_MoorePenrose_NoCent))
#' 
#' for (m in 1:5){
#'   cat("m = ", m, "\n")
#'   precision_higher_order_shrinkage_Cent = 
#'       ridge_higher_order_shrinkage(Y = t(X), m = m, centeredCov = TRUE, t = 10)
#'       
#'   precision_higher_order_shrinkage_NoCent = 
#'       ridge_higher_order_shrinkage(Y = t(X), m = m, centeredCov = FALSE, t = 10)
#'       
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_Cent$estimated_precision_matrix))
#'   
#'   print(FrobeniusLoss2(precision_higher_order_shrinkage_NoCent$estimated_precision_matrix))
#' }
#' 
#' 
#' precision_higher_order_shrinkage_Cent = 
#'       ridge_higher_order_shrinkage(Y = t(X), m = 1, centeredCov = TRUE, t = 100)
#' 
#' precision_target_identity_semioptimal_Cent = 
#'       ridge_target_identity_semioptimal(Y = t(X), centeredCov = TRUE, t = 100)
#'       
#' precision_target_general_semioptimal_Cent = 
#'       ridge_target_general_semioptimal(Y = t(X), centeredCov = TRUE, t = 100, Pi0 = diag(p))
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
#' FrobeniusLoss2(precision_higher_order_shrinkage_Cent$estimated_precision_matrix)
#' FrobeniusLoss2(precision_target_identity_semioptimal_Cent$estimated_precision_matrix)
#' 
#' 
#' @export
#' 
ridge_higher_order_shrinkage <- function(Y, m, centeredCov, t)
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
  
  # Regularized sample covariance matrix (Tikhonov regularization)
  S_t <- S + t * diag(nrow = p)
  
  S_t_inverse <- solve(S_t)
  
  q1 <- tr(S) / p
  q2 <- tr(S %*% S) / p - c_n * q1^2
  
  estimatedM = compute_M_t(m = m, c_n = c_n, q1 = q1, q2 = q2,
                           S_t_inverse = S_t_inverse,
                           t = t)
  
  # TODO: compute all estimators for smaller m here using submatrices of this matrix
  
  alpha = solve(estimatedM$M) %*% estimatedM$hm
  
  result = alpha[1] * Ip
  power_S_t_inverse = Ip
  
  for (k in 1:m){
    power_S_t_inverse = power_S_t_inverse %*% S_t_inverse
    result = result + alpha[k + 1] * power_S_t_inverse
  }
  
  return(list(
    estimated_precision_matrix = result,
    M = estimatedM$M,
    hm = estimatedM$hm,
    alpha = alpha,
    v = estimatedM$v
  ) )
}

