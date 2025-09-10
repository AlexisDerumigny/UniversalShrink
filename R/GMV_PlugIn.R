

#' Compute the optimal portfolio weights given an (estimated) precision matrix
#' 
#' This functions return the optimal portfolio weights defined by
#' \deqn{
#' w = \dfrac{\Sigma^{-1} \mathbf{1}}{\mathbf{1}^T \Sigma^{-1} \mathbf{1}}
#' }
#' i.e. \eqn{w = (w_1, \dots, w_p)} where for \eqn{j = 1, \dots, p}, we have
#' \deqn{
#' w_j = \dfrac{\sum_{k=1}^p [\Sigma^{-1}]_{j,k}}{
#' \sum_{k,l = 1}^p [\Sigma^{-1}]_{k,l}}.
#' }
#' Here \eqn{\mathbf{1}} is a vector \eqn{1} of size \eqn{p}, and
#' \eqn{\Sigma^{-1}} is the precision matrix (or an estimator thereof), of size
#' \eqn{p \times p}, where \eqn{p} is the number of assets that are considered
#' and \eqn{\Sigma} is the covariance matrix of the vector of asset returns.
#' 
#' 
#' @param estimatedPrecisionMatrix a matrix of size \eqn{p \times p}, which
#' is estimating the precision matrix \eqn{\Sigma^{-1}} 
#' ( = the inverse of the covariance matrix of the vector of asset returns that
#' are considered).
#' 
#' @returns a vector of size \eqn{p} with the computed portfolio weights, with
#' the minimum variance, i.e. the portfolio that has the lowest variance among
#' all portfolios.
#' 
#' @export
GMV_PlugIn <- function(estimatedPrecisionMatrix){
  
  # result = (estimatedPrecisionMatrix %*% ones) / 
  #                 (ones %*% estimatedPrecisionMatrix %*% ones)
  # Faster and equivalent expression:
  result = rowSums(estimatedPrecisionMatrix) / sum(estimatedPrecisionMatrix)
  
  return(result)
}

