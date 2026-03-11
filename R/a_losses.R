#' Quadratic losses of the estimator of a matrix
#' @name quadratic_losses
NULL


#' @export
#' @rdname quadratic_losses
NormFrobenius2 <- function(M, normalized){
  FrobNorm2 = tr( M %*% t(M) )
  if (normalized){
    p = ncol(M)
    return (FrobNorm2 / p)
  } else {
    return (FrobNorm2)
  }
}


#' @export
#' @rdname quadratic_losses
DistanceFrobenius2 <- function(M1, M2, normalized){
  FrobNorm2 = NormFrobenius2(M1 - M2, normalized = normalized)
  return (FrobNorm2)
}


#' This is a generic function to compute the Frobenius norm/loss of
#' (the estimator of) a matrix.
#' For a generic matrix \eqn{M}, the (squared) Frobenius norm is defined
#' in the following way
#' \deqn{
#' \code{NormFrobenius2}(\mathbf{M})=||\mathbf{M}||^2_F
#' = \text{\rm tr} \left[\mathbf{M}\mathbf{M}^\top\right]\,,
#' }
#' while (squared) Frobenius loss analogously is defined by
#' \deqn{
#' \code{LossFrobenius2}(\mathbf{M}, g(\boldsymbol{\Sigma}))
#' = ||\mathbf{M}-g(\mathbf{\boldsymbol{\Sigma}})||^2_F\,,
#' }
#' where \eqn{M} here denotes a suitable estimator and the function \eqn{g}
#' is \eqn{g(x)=x} when \code{type="covariance matrix"} is
#' chosen, otherwise \eqn{g(x)=1/x} for \code{type="precision matrix"}.
#' In the case \code{normalized=TRUE} (default), the above losses are
#' normalized by the matrix dimension \eqn{p}.
#' A related function for the case of precision matrices is
#' \code{LossInverseFrobenius2} which returns
#' \deqn{
#' \code{LossInverseFrobenius2}(\mathbf{M}, \Sigma)
#' = ||\mathbf{M} \boldsymbol{\Sigma} - I_p||^2_F\,,
#' }
#' where \eqn{I_p} is the identity matrix.
#' Furthermore we present alternative loss measures focused directly
#' on the spectra of the matrices: \code{DistanceEuclideanEigenvalues2} and 
#' \code{LossEuclideanEigenvalues2} defined as follows
#' \deqn{
#'  \code{DistanceEuclideanEigenvalues2}(\mathbf{M}_1, \mathbf{M}_2)=
#'  \sum\limits_{i=1}^{p} (\lambda_i(\mathbf{M}_1)-\lambda_i(\mathbf{M}_2))^2\,
#'  } and
#'  \deqn{
#'  \code{LossEuclideanEigenvalues2}(\mathbf{M}, g(\boldsymbol{\Sigma}) )= 
#'  \sum\limits_{i=1}^{p} 
#'  ( \lambda_i(\mathbf{M})-\lambda_i(g(\boldsymbol{\Sigma})) )^2
#' }
#' where \eqn{\lambda_i(\mathbf{A})} are the eigenvalues of a generic matrix
#' \eqn{\mathbf{A}}. Similarly, \code{normalized=TRUE} (by default) normalizes 
#' the losses, and the function \eqn{g(x)} is either \eqn{x} or \eqn{1/x},
#' which can be specified in  \code{type}. 
#'
#'
#' @param x,M An (estimated) square matrix of size \code{p}
#' @param M1,M2 two square matrices of the same dimension
#' 
#' @param Sigma,SigmaInv the true covariance matrix and its inverse
#' 
#' @param type target of the estimator. Must be a character vector of length 1
#' with one of the following:
#' \itemize{
#'   \item \code{type = "precision matrix"} corresponds to the
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{NormFrobenius2(x - SigmaInv)}, potentially normalized.
#'   
#'   \item \code{type = "covariance matrix"} corresponds to the
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{NormFrobenius2(x - Sigma )}, potentially normalized.
#' }
#' 
#' @param normalized if \code{TRUE}, the norm (or loss) is divided by the
#' matrix size \code{p}.
#' 
#' @param ... Additional arguments passed to methods.
#' 
#' 
#' @returns All these functions return a positive numeric value.
#' 
#' @examples
#' M = diag(c(1,2,3))
#' NormFrobenius2(M, normalized = TRUE)
#' NormFrobenius2(M, normalized = FALSE)
#' 
#' X <- MASS::mvrnorm(n = 100, mu = rep(0,3), Sigma = M)
#' estimatedCov_sample = cov(X)
#' 
#' # Losses for estimation of the covariance
#' LossFrobenius2(estimatedCov_sample, M, type = "covariance")
#' LossEuclideanEigenvalues2(estimatedCov_sample, M, type = "covariance")
#' 
#' 
#' 
#' @export
#' @rdname quadratic_losses
LossFrobenius2 <- function(x, Sigma, type, normalized = TRUE, ...) {
  UseMethod("LossFrobenius2")
}


#' @export
#' @rdname quadratic_losses
LossFrobenius2.matrix <- function(
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
      
      result = NormFrobenius2(x - SigmaInv, normalized = normalized)
    },
    
    "covariance matrix" = {
      result = NormFrobenius2(x - Sigma, normalized = normalized)
    },
    
    # default
    {
      stop("Type " , type, "is not implemented yet.")
    }
  )
  
  
  return (result)
}


#' @export
#' @rdname quadratic_losses
LossFrobenius2.EstimatedPrecisionMatrix <- function(
    x, Sigma, type = "precision matrix", normalized = TRUE, SigmaInv = NULL, ...)
{
  type = match.arg(type)
  
  if (type != "precision matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedPrecisionMatrix'.")
  }
  result = LossFrobenius2(x = as.matrix(x), Sigma = Sigma, SigmaInv = SigmaInv,
                          type = "precision matrix", normalized = normalized)
  
  return (result)
}



#' @export
#' @rdname quadratic_losses
LossFrobenius2.EstimatedCovarianceMatrix <- function(
    x, Sigma, type = "covariance matrix", normalized = TRUE, ...)
{
  type = match.arg(type)
  
  if (type != "covariance matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedCovarianceMatrix'.")
  }
  result = LossFrobenius2(as.matrix(x), Sigma = Sigma,
                          type = "covariance matrix", normalized = normalized)
  
  return (result)
}


#' @export
#' @rdname quadratic_losses
LossInverseFrobenius2 <- function(x, Sigma, normalized = TRUE, ...) {
  UseMethod("LossInverseFrobenius2")
}


#' @export
#' @rdname quadratic_losses
LossInverseFrobenius2.matrix <- function(x,
                                         Sigma,
                                         normalized = TRUE, ...)
{
  if (ncol(x) != nrow(x) || ncol(x) != nrow(Sigma) || ncol(x) != ncol(Sigma)){
    stop("x and Sigma should be square matrices of the same dimension. ",
         "Here dim(x) = c(", paste(dim(x), collapse = ","),
         ") and dim(Sigma) = c(", paste(dim(Sigma), collapse = ","), ")." )
  }
  p = ncol(Sigma)
  Ip = diag(p)
  
  result = NormFrobenius2(x %*% Sigma - Ip, normalized = normalized)
  
  return (result)
}


#' @export
#' @rdname quadratic_losses
LossInverseFrobenius2.EstimatedPrecisionMatrix <- function(
    x, Sigma, normalized = TRUE, ...)
{
  result = LossInverseFrobenius2(as.matrix(x), Sigma, normalized = normalized)
  
  return (result)
}


#' @export
#' @rdname quadratic_losses
DistanceEuclideanEigenvalues2 <- function(M1, M2, normalized){
  
  if (normalized){
    result = mean((eigen(M1)$values - eigen(M2)$values)^2)
  } else {
    result = sum((eigen(M1)$values - eigen(M2)$values)^2)
  }
  
  return (result)
}

#' @export
#' @rdname quadratic_losses
LossEuclideanEigenvalues2 <- function(x, Sigma, type,
                                      normalized = TRUE, ...) {
  UseMethod("LossEuclideanEigenvalues2")
}


#' @export
#' @rdname quadratic_losses
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
#' @rdname quadratic_losses
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
#' @rdname quadratic_losses
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




