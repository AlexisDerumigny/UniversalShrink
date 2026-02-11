

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
#' \doi{10.5705/ss.2012.328}
#' 
#' @examples
#' p = 200
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
#' precision_OptimalRidge = ridge_shrinkage_rescaled(t(X))
#' 
#' precisionTrue = solve(Sigma)
#' 
#' estimatedCov_NLshrink = cov_analytical_NL_shrinkage(X)
#' estimatedCov_QISshrink = cov_quadratic_inverse_shrinkage(X)
#' 
#' precision_NLshrink = solve(estimatedCov_NLshrink)
#' precision_QISshrink = solve(estimatedCov_QISshrink)
#' 
#' FrobeniusLoss2(precision_OptimalRidge, Sigma = Sigma)
#' FrobeniusLoss2(precision_NLshrink, Sigma = Sigma, type = "precision matrix")
#' FrobeniusLoss2(precision_QISshrink, Sigma = Sigma, type = "precision matrix")
#' 
#' 
#' @export
ridge_shrinkage_rescaled <- function (Y, eps = 1e-6, upp = pi/2 - 1e-6)
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
  
  S <- stats::cov(t(Y))
  
  p = nrow(Y)
  n = ncol(Y)
  c_n = p / n
  
  # Identity matrix of size p
  Ip = diag(nrow = p)
  
  # We are looking for an optimal rescaled ridge of the form
  # alpha * (Sn + t * Ip)^1
  # alpha can be optimized explicitly, but the objective has to be optimized
  # explicitly in t.
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
  
  hL_WPTZ_max <- stats::optim(1.5, hL_WPTZ, lower = eps, upper = upp, method = "L-BFGS-B",
                              control = list(fnscale = -1))
  t_optimal <- tan(hL_WPTZ_max$par)
  
  iS_t<-solve(S / t_optimal + Ip)
  tr_iS_t<-sum(diag(iS_t))/p
  a1<-1-tr_iS_t
  a2<-tr_iS_t-sum(diag(iS_t%*%iS_t))/p
  
  hR1<-a1/(1-c_n*a1)
  hR2<-a1/((1-c_n*a1)^3)-a2/((1-c_n*a1)^4)
  alpha_optimal <- hR1 / hR2
  
  # We are looking for an optimal rescaled ridge of the form
  # alpha * (Sn + t * Ip)^1
  # This is the final estimator of the precision matrix
  # using the optimal parameters that were found.
  iS_WPTZ <- alpha_optimal * solve(S + t_optimal * Ip)
  
  result = list(
    estimated_precision_matrix = iS_WPTZ
  )
  
  class(result) <- c("EstimatedPrecisionMatrix")
  
  return (result)
}

