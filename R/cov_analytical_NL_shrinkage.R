

#' Analytical non-linear shrinkage from Ledoit and Wolf (2020)
#' 
#' This function estimates the covariance matrix as a given data set using 
#' classical nonlinear shrinkage from Ledoit and Wolf (2020).
#' 
#' For \eqn{i\in\{1,\ldots,p\}}, the nonlinear shrinkage estimator of 
#' Ledoit and Wolf is given by
#' \deqn{
#' \mathbf{S}_{NLSh}= \mathbf{U} \text{diag}(d_1^{or},...,d_p^{or})
#' \mathbf{U}^\top, \quad d_i^{or}=\left\{
#'  \begin{array}{ll}
#'  \frac{d_i}{|1-c-c d_i \breve{m}_{F}(d_i)|^2},    & \text{if}~~ d_i>0,\\
#'  \frac{1}{(c-1)\breve{m}_{\underline{F}}(0)}, & \text{if}~~ d_i=0,
#'  \end{array}
#'  \right.
#'  }
#'  where \eqn{\mathbf{U}=(\mathbf{u}_1,...\mathbf{u}_p)} is the matrix with the
#'   sample eigenvectors of \eqn{\mathbf{S}_n}, \eqn{d_i}, \eqn{i=1,\ldots,} are
#'   the sample eigenvalues of \eqn{\mathbf{S}_n} and 
#'   \eqn{\breve{m}_{F}(x)=\lim\limits_{z\to x}m_{F}(z)} with \eqn{m_{F}(z)} the
#'    limiting Stieltjes transform of the sample covariance matrix. A numerical 
#'    approach to estimate \eqn{\breve{m}_{F}(x)} is provided in Ledoit and Wolf
#'     (2020) and is available in the R-package \eqn{\textit{HDShOP}}. 
#' 
#' 
#' @param X data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @returns the estimator of the covariance matrix
#' (a `p` by `p` matrix).
#' TODO: update this
#' 
#' @references 
#' Ledoit, O., & Wolf, M. (2020).
#' Analytical nonlinear shrinkage of large-dimensional covariance matrices.
#' The Annals of Statistics, 48(5), 3043-3065.
#' \doi{10.1214/19-AOS1921}
#' 
#' Bodnar, T., S. Dmytriv, Y. Okhrin, D. Otryakhin, & N. Parolya. 2024. HDShOP: 
#' High-Dimensional Shrinkage Optimal Portfolios. 
#' R package version 0.1.7. \doi{10.32614/CRAN.package.HDShOP}
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
#' FrobeniusLoss2(estimatedCov_sample, Sigma, type = "covariance")
#' FrobeniusLoss2(estimatedCov_shrink, Sigma)
#' 
#' LossEuclideanEigenvalues2(estimatedCov_sample, Sigma, type = "covariance")
#' LossEuclideanEigenvalues2(estimatedCov_shrink, Sigma)
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
  
  estimatedSigma = u %*% diag(dtilde) %*% t(u)
  
  result = list(
    estimated_covariance_matrix = estimatedSigma
  )
  
  class(result) <- c("EstimatedCovarianceMatrix")
  
  return (result)
}



