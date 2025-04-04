

#' Optimal ridge shrinkage for estimation of the precision matrix
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param eps,upp search interval for the best penalization parameter
#' (used in the numerical optimization).
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix, i.e. the inverse of the covariance matrix).
#' 
#' @references 
#' Wang, C., Pan, G., Tong, T., & Zhu, L. (2015).
#' Shrinkage estimation of large dimensional precision matrix using random matrix theory.
#' Statistica Sinica, 993-1008.
#' \link{http://dx.doi.org/10.5705/ss.2012.328}
#' 
#' @examples
#' p = 200
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
#' precision_OptimalRidge = optimal_ridge_shrinkage(t(X))
#' 
#' precisionTrue = solve(Sigma)
#' 
#' estimatedCov_NLshrink = analytical_NL_shrinkage(t(X))
#' estimatedCov_QISshrink = quadratic_inverse_shrinkage(X)
#' 
#' precision_NLshrink = solve(estimatedCov_NLshrink)
#' precision_QISshrink = solve(estimatedCov_QISshrink)
#' 
#' mean((eigen(precisionTrue)$values - eigen(precision_OptimalRidge)$values)^2)
#' mean((eigen(precisionTrue)$values - eigen(precision_NLshrink)$values)^2)
#' mean((eigen(precisionTrue)$values - eigen(precision_QISshrink)$values)^2)
#' 
#' 
#' @export
optimal_ridge_shrinkage <- function (Y, eps = 1e-6, upp = pi/2 - 1e-6)
{
  if (eps <= 0 || eps > pi/2){
    stop("'eps' must be between 0 and pi/2")
  }
  if (upp <= 0 || upp > pi/2){
    stop("'upp' must be between 0 and pi/2")
  }
  if (eps >= upp){
    stop("'eps' must be strictly smaller than 'upp'.")
  }
  
  S <- cov(t(Y))
  
  p = nrow(Y)
  n = ncol(Y)
  c_n = p / n
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  hL_WPTZ<- function(u)
  {
    t<-tan(u)
    iS_t<-solve(S/t + Ip)
    tr_iS_t<-sum(diag(iS_t))/p
    a1<-1-tr_iS_t
    a2<-tr_iS_t-sum(diag(iS_t%*%iS_t))/p
    
    hR1<-a1/(1-c_n*a1)
    hR2<-a1/((1-c_n*a1)^3)-a2/((1-c_n*a1)^4)
    
    hL_WPTZ<-hR1^2/hR2
    return(hL_WPTZ)
  }
  
  hL_WPTZ_max <- optim(1.5, hL_WPTZ, lower = eps, upper = upp, method = "L-BFGS-B",
                       control = list(fnscale = -1))
  hL_WPTZ_bet <- tan(hL_WPTZ_max$par)
  
  iS_t<-solve(S/hL_WPTZ_bet + Ip)
  tr_iS_t<-sum(diag(iS_t))/p
  a1<-1-tr_iS_t
  a2<-tr_iS_t-sum(diag(iS_t%*%iS_t))/p
  
  hR1<-a1/(1-c_n*a1)
  hR2<-a1/((1-c_n*a1)^3)-a2/((1-c_n*a1)^4)
  hL_WPTZ_alp<-hR1/hR2
  
  iS_WPTZ <- hL_WPTZ_alp*solve(S +hL_WPTZ_bet * Ip)
  
  return (iS_WPTZ)
}

