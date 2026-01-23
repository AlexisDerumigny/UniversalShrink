

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
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @export
GMV_Moore_Penrose <- function(Y, centeredCov = TRUE)
{
  iS_MP = Moore_Penrose(Y = Y, centeredCov = centeredCov)
  GMV_MP = GMV_PlugIn(iS_MP)
  
  return (GMV_MP)
}


