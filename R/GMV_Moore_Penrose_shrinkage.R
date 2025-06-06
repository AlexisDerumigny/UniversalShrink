

#' First-order shrinkage of the Moore-Penrose portfolio towards a general target
#'
#'
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
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
#'   GMV_Moore_Penrose_shrinkage(Y = t(X), centeredCov = TRUE)
#' 
#' outOfSampleVariance = t(GMV_MP_shrinkage_Cent) %*% Sigma %*% GMV_MP_shrinkage_Cent
#' 
#' ones = rep(1, length = p)
#' V_GMV = 1 / ( t(ones) %*% solve(Sigma) %*% ones)
#' 
#' Loss_GMV_Moore_Penrose_shrinkage = (outOfSampleVariance - V_GMV) / V_GMV
#' 
#' GMV_MP_Cent = GMV_Moore_Penrose(Y = t(X), centeredCov = TRUE)
#' outOfSampleVariance = t(GMV_MP_Cent) %*% Sigma %*% GMV_MP_Cent
#' 
#' Loss_GMV_Moore_Penrose = (outOfSampleVariance - V_GMV) / V_GMV
#' 
#' # Shrinkage helps to reduce the loss
#' stopifnot(Loss_GMV_Moore_Penrose_shrinkage < Loss_GMV_Moore_Penrose)
#' 
#' 
#' 
#' @export
GMV_Moore_Penrose_shrinkage <- function(Y, b = NULL, centeredCov = TRUE){
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
  
  bip<-matrix(rep(1,p),p,1)
  tbip<-t(bip)
  
  bipiSbip<-sum(tbip%*%iS_MP%*%bip)/p
  bipiS2bip<-sum(tbip%*%iS_MP%*%iS_MP%*%bip)/p
  bipiS3bip<-sum(tbip%*%iS_MP%*%iS_MP%*%iS_MP%*%bip)/p
  
  bipSbip<-sum(tbip%*%S%*%bip)/p
  
  hv0 <- c_n * trS1
  
  d0<-1-tbip%*%S%*%iS_MP%*%bip/p
  d1<-bipiSbip/p/trS2/c_n  
  d1_b<-((1-d0)/(hv0)-d1)/hv0
  d3<-(bipiS3bip/trS2^3+2*bipiSbip*trS3^2/p/trS2^5-(bipiS2bip+trS4*bipiSbip)/trS2^4)/c_n^3
  
  alp_ShMP <- sum(bipSbip - d1_b / d1) / sum(bipSbip-2*d1_b/d1+d3/d1^2)
  
  cat("alp_ShMP = ", alp_ShMP, "\n")
  
  w_ShMP <- alp_ShMP*w_MP+(1-alp_ShMP)*bip/p
  
  return (w_ShMP)
}
