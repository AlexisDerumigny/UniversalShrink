

#' Optimal portfolio weights based on a plug-in of the Moore-Penrose inverse
#' of the sample covariance matrix
#' 
#' For an \eqn{n \times p} observation matrix \eqn{\mathbf{X}_n} and the
#' corresponding sample covariance matrix \eqn{\mathbf{S}_n} of size
#' \eqn{p \times p}, this estimator is given by
#' \deqn{
#' \mathbf{w}_{MP}=\dfrac{\mathbf{S}^+_n\mathbf{1}}{\mathbf{1}^\top
#' \mathbf{S}^+_n\mathbf{1}},
#' }
#' where \eqn{\mathbf{S}^+_n} is the Moore-Penrose (pseudo) inverse of 
#' \eqn{\mathbf{S}_n} and \eqn{\mathbf{1}} is a \eqn{p}-dimensional vector of
#' ones.
#' 
#'
#' @param X data matrix (rows are observations, columns are features).
#' @inheritParams Moore_Penrose
#' 
#' @returns a vector of size \eqn{p} giving the optimal portfolio weights, with
#' the elements summing up to \eqn{1}.
#' 
#' @seealso \code{\link{LossRelativeOutOfSampleVariance}} for computing the loss
#' of the portfolio.
#' 
#' @examples
#' Sigma = diag(1:5)
#' 
#' X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = Sigma)
#' weights100 = GMV_Moore_Penrose(X)
#' 
#' X <- MASS::mvrnorm(n = 3, mu = rep(0,5), Sigma = Sigma)
#' weights3 = GMV_Moore_Penrose(X)
#' 
#' LossRelativeOutOfSampleVariance(weights100, Sigma)
#' LossRelativeOutOfSampleVariance(weights3, Sigma)
#' 
#' trueWeights = rowSums(solve(Sigma)) / sum(solve(Sigma))
#' 
#' NormFrobenius2(trueWeights - weights100, normalized = FALSE)
#' NormFrobenius2(trueWeights - weights3, normalized = FALSE)
#' 
#' @export
GMV_Moore_Penrose <- function(X, centeredCov = TRUE)
{
  iS_MP = Moore_Penrose(X = X, centeredCov = centeredCov)
  GMV_MP = GMV_PlugIn(iS_MP)
  
  return (GMV_MP)
}


