

# Trace function of a matrix
tr <- function(M){
  if (is.matrix(M)){
    result = sum(diag(M))
  } else {
    n = nrow(M)
    diagM = M[cbind(1:n, 1:n)]
    result = sum(diagM)
  }
  return (result)
}


#' Conversion of estimated matrices to matrix class
#' 
#' @param x object to be converted
#' @param ... other arguments passed from methods, currently ignored.
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


#' Constructor for warning conditions of the package
#'
#' @noRd
UniversalShrink_warning_condition_base <- function(message, subclass = NULL,
                                                   call = sys.call(-1), ...) {
  # warningCondition() automatically adds 'warning' and 'condition' to the class
  return (
    warningCondition(
      message = message,
      class = c(subclass, "UniversalShrinkWarning"), # We add a base warning class
      call = call,
      ... # Allows for additional custom fields
    )
  )
}
