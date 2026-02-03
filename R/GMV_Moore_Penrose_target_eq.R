

GMV_Moore_Penrose_target_eq <- function(Y, centeredCov = TRUE, verbose = 2){
  
  # Get sizes of Y
  p = nrow(Y)
  n = ncol(Y)
  c_n = concentr_ratio(n = n, p = p, centeredCov = centeredCov, verbose = verbose)
  
  # Sample covariance matrix
  S <- cov_with_centering(X = t(Y), centeredCov = centeredCov)
  
  # Moore-Penrose inverse of the sample covariance matrix
  iS_MP <- as.matrix(Moore_Penrose(Y = Y, centeredCov = centeredCov))
  
  w_MP = GMV_PlugIn(estimatedPrecisionMatrix = iS_MP)
  
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
  
  return (w_ShMP)
}
