

# Trace function of a matrix
tr <- function(M){
  return (sum(diag(M)))
}


# Compute the matrix M for the higher-order shrinkage
compute_M <- function(m, n, p, ihv0, D_MP, q1 = q1, q2 = q2)
{
  v <- rep(0, 2*m)
  h <- rep(0, 2*m+1)
  d <- matrix(0, nrow = 2 * m, ncol = 2)
  s <- rep(q2 , (2*m+1))
  
  h[2] <- h2
  h[3] <- h3
  
  for (i in 1:(2*m))
  {
    v[i] <- (-1)^i * factorial(i) * c_n * (1/p) * tr(D_MP^(i+1))
  }
  
  
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
        d[k,l] <- ihv0 * (ihv0 * (q1 - ihv0 * (1 / c_n)) - d[1,1])
      } else if (l == 2 && k > 1) {
        d[k,l] <- ihv0 * (d[k-1,2] - d[k,1])
      }
    }
  }
  
  for (j in 1:(2*m))
  {
    ss <- 0
    for (k in 1:j)
    {
      ss <- ss + (-1)^(j+k+1) * factorial(k) / factorial(j) * d[k,2] * 
        kStatistics::e_eBellPol(j, k, c(v[1:(j-k+1)], rep(0,k-1)))
    }
    s[j+1] <- ss
  }
  
  M <- s[1:(m+1)]
  for (j in 2:(m+1))
  {
    M <- cbind(M, s[j:(m+j)])
  }
  
  M <- apply(M, c(1,2), Re)
  
  return(M)
}



#' Higher order shrinkage
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param m order of the shrinkage. Should be at least 1.
#' 
#' 
higher_order_shrinkage <- function(Y, m, centeredCov)
{
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = p / n
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # Sample covariance matrix
    S <- Y %*% Jn %*% t(Y) / n
    
    # We remove the last eigenvector because the eigenvalues are sorted
    # in decreasing order.
    Hn = eigen(Jn)$vectors[, -n]
    Ytilde = Y %*% Hn
    
    # Inverse companion covariance
    iYtilde <- solve(t(Ytilde) %*% Ytilde / n)
    
    # Moore-Penrose inverse
    iS_MP <- Ytilde %*% iYtilde %*% iYtilde %*% t(Ytilde)/n
    
  } else {
    S <- Y %*% t(Y)/n
    
    # Inverse companion covariance
    iY <- solve(t(Y) %*% Y / n)
    
    # Moore-Penrose inverse
    iS_MP <- Y %*% iY %*% iY %*% t(Y)/n
  }
  
  D_MP <- diag(eigen(iS_MP)$values)
  
  trS1 <- tr(iS_MP) / p
  trS2 <- tr(iS_MP %*% iS_MP) / p
  trS3 <- tr(iS_MP %*% iS_MP %*% iS_MP) / p
  
  hv0 <- c_n * trS1
  ihv0 <- 1 / hv0
  ihv0_2 <- ihv0^2
  ihv0_3 <- ihv0^3
  
  ##### Shrinkage m order
  
  
  hm <- function(m)
  {
    hm <- rep(0, m + 1)
    h  <- rep(0, m + 1)
    d  <- rep(0, m)
    v  <- rep(0, m)
    
    
    for (i in 1:m) {
      v[i]<- (-1)^i * factorial(i) * c_n * (1/p) * tr(D_MP^(i+1))
    }
    
    hm[1]<- q1
    
    h[2]<- h2
    h[3]<- h3
    
    if (m > 1) { 
      for (i in 3:m) {
        sh <- 0
        for (k in 2:(i-1)) {
          sh <- sh + (-1)^k * factorial(k) * h[k+1] * 
            kStatistics::e_eBellPol(i, k, c(v[1:(i-k+1)], rep(0,k-1)))
        }
        h[i+1] <- (v[i] + v[1] * sh) / ((v[1])^(i+1) * (-1)^(i+1) * factorial(i))
      }
    }
    
    for (k in 1:m)
    {
      d[k] <- 1 / c_n * (1 / (hv0^(k+1)) - h[k+1])
    }
    
    hm[1]<-q1
    for (j in 1:m)
    {
      s<-0
      for (k in 1:j)
      {
        s<- s + (-1)^(j+k+1)*factorial(k)/factorial(j)*d[k] *
          kStatistics::e_eBellPol(j,k,c(v[1:(j-k+1)],rep(0,k-1)))
      }
      hm[j+1]<- s
    }
    hm<- Re(hm)
    return(hm)
  }
  
  
  # alpha_m1 <- solve(compute_M(m = 1, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
  #                             q1 = q1, q2 = q2)) %*% hm(1)
  
  alpha_m2 <- solve(compute_M(m = 2, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
                              q1 = q1, q2 = q2)) %*% hm(2)
  
  alpha_m3 <- solve(compute_M(m = 3, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
                              q1 = q1, q2 = q2)) %*% hm(3)
  
  alpha_m4 <- solve(compute_M(m = 4, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
                              q1 = q1, q2 = q2)) %*% hm(4)
  
  alpha_m5 <- solve(compute_M(m = 5, n = n, p = p, ihv0 = ihv0, D_MP = D_MP,
                              q1 = q1, q2 = q2)) %*% hm(5)
  
  
  # high_shrink1<- alpha_m1[1]*Ip+alpha_m1[2]*iS_MP
  high_shrink2 <- alpha_m2[1]*Ip + alpha_m2[2]*iS_MP + alpha_m2[3]*iS_MP%*%iS_MP
  
  high_shrink3 <- alpha_m3[1]*Ip + alpha_m3[2]*iS_MP + alpha_m3[3]*iS_MP%*%iS_MP + 
    alpha_m3[4]*iS_MP%*%iS_MP%*%iS_MP
  
  high_shrink4 <- alpha_m4[1]*Ip + alpha_m4[2]*iS_MP + alpha_m4[3]*iS_MP%*%iS_MP + 
    alpha_m4[4]*iS_MP%*%iS_MP%*%iS_MP + alpha_m4[5]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP
  
  high_shrink5 <- alpha_m5[1]*Ip + alpha_m5[2]*iS_MP + alpha_m5[3]*iS_MP%*%iS_MP + 
    alpha_m5[4]*iS_MP%*%iS_MP%*%iS_MP + alpha_m5[5]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP + 
    alpha_m5[6]*iS_MP%*%iS_MP%*%iS_MP%*%iS_MP%*%iS_MP
  
  
}

