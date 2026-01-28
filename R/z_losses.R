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
#' @param x,M An (estimated) square matrix of size \code{p}
#' @param M1,M2 two square matrices of the same dimension
#' 
#' @param Sigma,SigmaInv the true covariance matrix and its inverse
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
#' @param portfolioWeights the vector of weights of a given portfolio, of which
#' we want to determine the loss
#' 
#' @param normalized if \code{TRUE}, the Frobenius norm is divided by the matrix
#' size \code{p}.
#' 
#' @param ... Additional arguments passed to methods.
#' 
#' 
#' @returns All these functions return a positive numeric value.
#' \code{FrobeniusNorm2} returns the squared Frobenius norm of a matrix.
#' \code{FrobeniusLoss2} returns the (normalized) Frobenius loss.
#'
#' \code{LossRelativeOutOfSampleVariance} returns a positive numeric value,
#' with attributes \code{"V_portfolio"} and \code{"V_GMV"}, which are respectively
#' the (out of sample) variances of the given \code{portfolioWeights} and of the
#' GMV portfolio.
#' 
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
  
  if (missing(type)){
    stop("'type' must be specified. Either 'precision' or 'covariance'.")
  }
  
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
  result = FrobeniusLoss2(as.matrix(x), Sigma,
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
  result = FrobeniusLoss2(as.matrix(x), Sigma,
                          type = "covariance matrix", normalized = normalized)
  
  return (result)
}


#' @export
#' @rdname FrobeniusLoss2
DistanceEuclideanEigenvalues2 <- function(M1, M2, normalized){
  
  if (normalized){
    result = mean((eigen(M1)$values - eigen(M2)$values)^2)
  } else {
    result = sum((eigen(M1)$values - eigen(M2)$values)^2)
  }
  
  return (result)
}

#' @export
#' @rdname FrobeniusLoss2
LossEuclideanEigenvalues2 <- function(x, Sigma, type,
                                      normalized = TRUE, ...) {
  UseMethod("LossEuclideanEigenvalues2")
}


#' @export
#' @rdname FrobeniusLoss2
LossEuclideanEigenvalues2.matrix <- function(
    x,
    Sigma,
    type = c("precision matrix", "covariance matrix"),
    normalized = TRUE, 
    SigmaInv = NULL, ...)
{
  if (ncol(x) != nrow(x) || ncol(x) != nrow(Sigma) || ncol(x) != ncol(Sigma)){
    stop("x and Sigma should be square matrices of the same dimension. ",
         "Here dim(x) = c(", paste(dim(x), collapse = ","),
         ") and dim(Sigma) = c(", paste(dim(Sigma), collapse = ","), ")." )
  }
  p = ncol(Sigma)
  
  if (missing(type)){
    stop("'type' must be specified. Either 'precision' or 'covariance'.")
  }
  
  type = match.arg(type)
  
  switch (
    type,
    
    "precision matrix" = {
      if (is.null(SigmaInv)) {
        SigmaInv = solve(Sigma)
      }
      result = DistanceEuclideanEigenvalues2(x, SigmaInv, normalized = normalized)
    },
    
    "covariance matrix" = {
      result = DistanceEuclideanEigenvalues2(x, Sigma, normalized = normalized)
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
LossEuclideanEigenvalues2.EstimatedPrecisionMatrix <- function(
    x, Sigma, type = "precision matrix", normalized = TRUE, SigmaInv = NULL, ...)
{
  type = match.arg(type)
  
  if (type != "precision matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedPrecisionMatrix'.")
  }
  result = LossEuclideanEigenvalues2(
    as.matrix(x), Sigma = Sigma,
    type = "precision matrix", normalized = normalized, SigmaInv = SigmaInv)
  
  return (result)
}


#' @export
#' @rdname FrobeniusLoss2
LossEuclideanEigenvalues2.EstimatedCovarianceMatrix <- function(
    x, Sigma, type = "covariance matrix", normalized = TRUE, ...)
{
  type = match.arg(type)
  
  if (type != "covariance matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedCovarianceMatrix'.")
  }
  result = LossEuclideanEigenvalues2(
    as.matrix(x), Sigma,
    type = "covariance matrix", normalized = normalized)
  
  return (result)
}


#' @export
#' @rdname FrobeniusLoss2
LossRelativeOutOfSampleVariance <- function(portfolioWeights, Sigma, SigmaInv = NULL){
  if (is.null(SigmaInv)){
    SigmaInv = solve(Sigma)
  }
  p = length(portfolioWeights)
  
  outOfSampleVariance = t(portfolioWeights) %*% Sigma %*% portfolioWeights

  ones = rep(1, length = p)
  V_GMV = 1 / ( t(ones) %*% SigmaInv %*% ones)

  Loss = as.numeric( (outOfSampleVariance - V_GMV) / V_GMV )
  
  attr(Loss, "V_portfolio") <- outOfSampleVariance
  attr(Loss, "V_GMV") <- V_GMV
  
  return (Loss)
}

