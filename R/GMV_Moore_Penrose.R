

#' Optimal portfolio weights based on a plug-in of the Moore-Penrose inverse
#' of the sample covariance matrix
#' 
#' Having the obervation matrix \eqn{\mathbf{Y}_n} and the corresponding sample 
#' covariance matrix \eqn{\mathbf{S}_n} this estimator is given 
#' by
#' \deqn{
#' \mathbf{w}_{MP}=\frac{\mathbf{S}^+_n\mathbf{1}}{\mathbf{1}^\top\mathbf{S}^+_n\mathbf{1}},
#' }
#' where \eqn{\mathbf{S}^+_n} is the Moore-Penrose (pseudo) inverse of 
#'\eqn{\mathbf{S}_n} and and \eqn{\mathbf{1}} is a vector of ones.
#' 
#'
#' @param X data matrix (rows are observations, columns are features).
#' @inheritParams Moore_Penrose
#' 
#' @export
GMV_Moore_Penrose <- function(X, centeredCov = TRUE)
{
  iS_MP = Moore_Penrose(X = X, centeredCov = centeredCov)
  GMV_MP = GMV_PlugIn(iS_MP)
  
  return (GMV_MP)
}


