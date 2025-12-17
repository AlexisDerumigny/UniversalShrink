

#' Analytical non-linear shrinkage from Ledoit and Wolf (2020)
#' 
#' 
#' This function estimates the covariance matrix as a given data set using
#' non-linear shrinkage.
#' 
#' @param x data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @returns the estimator of the covariance matrix
#' (a `p` by `p` matrix).
#' 
#' @references 
#' Ledoit, O., & Wolf, M. (2020).
#' Analytical nonlinear shrinkage of large-dimensional covariance matrices.
#' The Annals of Statistics, 48(5), 3043-3065.
#' \link{https://doi.org/10.1214/19-AOS1921}
#' 
#' @examples
#' p = 200
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = 100, mu = mu, Sigma=Sigma)
#' estimatedCov_sample = cov(X)
#' estimatedCov_shrink = cov_analytical_NL_shrinkage(t(X))
#' 
#' # We now compare the distance between the true and both estimators.
#' mean((eigen(Sigma)$values - eigen(estimatedCov_sample)$values)^2)
#' mean((eigen(Sigma)$values - eigen(estimatedCov_shrink)$values)^2)
#' 
#' @export
cov_analytical_NL_shrinkage = function(x){
  # the original version suggested that p is # of columns
  p = nrow(x)
  n = ncol(x)
  sampleC = stats::cov(t(x))
  eig = eigen(sampleC)
  u = eig$vectors[,p:1]
  lambda = rev(eig$values)
  
  if(p <= n){
    lambda = lambda[max(1, p-n+1):p]
    L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
    h = n^(-1/3)
    H = h * t(L)
    x <- (L - t(L)) / H # This is a different x than before
    ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
    Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
    Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
    Hftilde = rowMeans(Hftemp / H)
    dtilde = lambda / ((pi*(p/n)*lambda*ftilde)^2 + (1-(p/n)-pi*(p/n)*lambda*Hftilde)^2);
  } else{
    lambda = lambda[max(1, p-n+2):p]
    L = matrix(rep(lambda, min(p, n-1)), nrow = length(lambda))
    h = n^(-1/3)
    H = h * t(L)
    x <- (L - t(L)) / H # This is a different x than before
    ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
    Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
    Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
    Hftilde = rowMeans(Hftemp / H)
    
    Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) * log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1/lambda)
    dtilde0 = 1/(pi*(p-n)/n * Hftilde0)
    dtilde1 = lambda/(pi^2*lambda^2 * (ftilde^2 + Hftilde^2))
    dtilde = c(dtilde0 * rep(1, p-n+1), dtilde1)
  }
  u %*% diag(dtilde) %*% t(u)
}



