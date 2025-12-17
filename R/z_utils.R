

# Trace function of a matrix
tr <- function(M){
  return (sum(diag(M)))
}



#' @export
#' @rdname FrobeniusLoss2
FrobeniusNorm2 <- function(M, normalized){
  FrobNorm2 = tr( M %*% t(M) )
  if (normalized){
    p = ncol(M)
    return (FrobNorm2 / p)
  } else {
    return (FrobNorm2)
  }
}


#' Frobenius norm and Frobenius losses of the estimator of a matrix
#'
#' Generic function to calculate the Frobenius norm/loss of (the estimator of) a
#' matrix.
#'
#' @param x,M An (estimated) square matrix of size \code{p},
#' 
#' @param Sigma the true covariance matrix
#' 
#' @param type target of the estimator. Must be a character vector of length 1
#' with one of the following:
#' \itemize{
#'   \item \code{type = "precision matrix"} corresponds to the normalized
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{FrobeniusNorm2(x \%*\% Sigma - diag(p) ) / p}.
#'   
#'   \item \code{type = "covariance matrix"} corresponds to the normalized
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{FrobeniusNorm2(x - Sigma ) / p}.
#' }
#' 
#' @param normalized if \code{TRUE}, the Frobenius norm is divided by the matrix
#' size \code{p}.
#' 
#' @param ... Additional arguments passed to methods.
#' 
#' @returns All these functions return a positive numeric value.
#' \code{FrobeniusNorm2} returns the squared Frobenius norm of a matrix.
#' \code{FrobeniusLoss2} returns the (normalized) Frobenius loss.
#'
#' @export
FrobeniusLoss2 <- function(x, Sigma, type, normalized = TRUE, ...) {
  UseMethod("FrobeniusLoss2")
}


#' @export
#' @rdname FrobeniusLoss2
FrobeniusLoss2.matrix <- function(x,
                                  Sigma,
                                  type = c("precision matrix", "covariance matrix"),
                                  normalized = TRUE, ...)
{
  if (ncol(x) != nrow(x) || ncol(x) != nrow(Sigma) || ncol(x) != ncol(Sigma)){
    stop("x and Sigma should be square matrices of the same dimension. ",
         "Here dim(x) = c(", paste(dim(x), collapse = ","),
         ") and dim(Sigma) = c(", paste(dim(Sigma), collapse = ","), ")." )
  }
  p = ncol(Sigma)
  
  type = match.arg(type)
  
  switch (
    type,
    
    "precision matrix" = {
      result = FrobeniusNorm2(x %*% Sigma - diag(p), normalized = normalized)
    },
    
    "covariance matrix" = {
      result = FrobeniusNorm2(x - Sigma, normalized = normalized)
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
FrobeniusLoss2.EstimatedPrecisionMatrix <- function(
    x, Sigma, type = "precision matrix", normalized = TRUE, ...)
{
  type = match.arg(type)
  
  if (type != "precision matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedPrecisionMatrix'.")
  }
  result = FrobeniusLoss2(x$estimated_precision_matrix, Sigma,
                          type = "precision matrix", normalized = normalized)
  
  return (result)
}



#' @export
#' @rdname FrobeniusLoss2
FrobeniusLoss2.EstimatedCovarianceMatrix <- function(
    x, Sigma, type = "covariance matrix", normalized = TRUE, ...)
{
  type = match.arg(type)
  
  if (type != "covariance matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedCovarianceMatrix'.")
  }
  result = FrobeniusLoss2(x$estimated_covariance_matrix, Sigma,
                          type = "covariance matrix", normalized = normalized)
  
  return (result)
}

