

#' Compute the optimal portfolio weights given an (estimated) precision matrix
#' 
#' @param estimatedPrecisionMatrix a matrix of size \code{p * p}, which
#' is estimating the precision matrix ( = the inverse of the covariance matrix)
#' 
#' @returns a vector of size \code{p} with the computed portfolio weights
#' 
#' @export
GMV_PlugIn <- function(estimatedPrecisionMatrix){
  
  # result = (estimatedPrecisionMatrix %*% ones) / (ones %*% estimatedPrecisionMatrix %*% ones)
  # Faster and equivalent expression:
  result = rowSums(estimatedPrecisionMatrix) / sum(estimatedPrecisionMatrix)
  
  return(result)
}

