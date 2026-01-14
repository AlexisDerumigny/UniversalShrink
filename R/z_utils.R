

# Trace function of a matrix
tr <- function(M){
  return (sum(diag(M)))
}


#' Conversion of estimated matrices to matrix class
#' 
#' @param x object to be converted
#' 
#' @return the underlying estimated matrix
#' 
#' @export
as.matrix.EstimatedPrecisionMatrix <- function(x, ...){
  return (x$estimated_precision_matrix)
}

#' @rdname as.matrix.EstimatedPrecisionMatrix
#' @export
as.matrix.EstimatedCovarianceMatrix <- function(x, ...){
  return (x$estimated_covariance_matrix)
}

