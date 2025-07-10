

# Trace function of a matrix
tr <- function(M){
  return (sum(diag(M)))
}


#' Squared Frobenius norm of a matrix
#' 
#' @export
FrobeniusNorm2 <- function(M){
  return ( tr( M %*% t(M) ) )
}


#' Frobenius loss of the estimator of a precision matrix
#'
#' Generic function to calculate the Frobenius loss of the estimator of a
#' precision matrix
#'
#' @param x An (estimated) matrix, either estimated precision matrix or estimated
#' covariance matrix
#' @param Sigma the true covariance matrix
#' @param type target of the estimator. Must be a character vector of length 1
#' with one of the following:
#' \itemize{
#'   \item \code{type = "precision matrix"} corresponds to the normalized
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{FrobeniusNorm2(x \%*\% Sigma - diag(p) ) / p}.
#' }
#' @param ... Additional arguments passed to methods.
#'
#' @export
FrobeniusLoss2 <- function(x, Sigma, type = "precision matrix", ...) {
  UseMethod("FrobeniusLoss2")
}


#' @export
#' @rdname FrobeniusLoss2
FrobeniusLoss2.matrix <- function(x, Sigma, type = "precision matrix", ...) {
  if (ncol(x) != nrow(x) || ncol(x) != nrow(Sigma) || ncol(x) != ncol(Sigma)){
    stop("x and Sigma should be square matrices of the same dimension. ",
         "Here dim(x) = c(", paste(dim(x), collapse = ","),
         ") and dim(Sigma) = c(", paste(dim(Sigma), collapse = ","), ")." )
  }
  p = ncol(Sigma)
  
  switch (
    type,
    
    "precision matrix" = {
      result = FrobeniusNorm2(x %*% Sigma - diag(p) ) / p
    },
    
    # default
    {
      stop("Type " , type, "is not implemented yet.")
    }
  )
  
  
  return (result)
}


#' @export
#' @rdname FrobeniusLoss2
FrobeniusLoss2.EstimatedPrecisionMatrix <- function(x, Sigma,
                                                    type = "precision matrix", ...) {
  if (type != "precision matrix"){
    stop("Type is chosen to be ", type, "but x is of type", EstimatedPrecisionMatrix)
  }
  result = FrobeniusLoss2(x$estimated_precision_matrix, Sigma, type = "precision matrix")
  
  return (result)
}


