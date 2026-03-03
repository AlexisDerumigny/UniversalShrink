#' Quadratic losses of the estimator of a matrix or portfolio weights
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


#' Generic function to calculate the Frobenius norm/loss of (the estimator of) a
#' matrix. For a generic matrix \eqn{M}, the (squared) Frobenius norm is defined in the following way
#' \deqn{
#' \code{NormFrobenius2}(\mathbf{M})=||\mathbf{M}||^2_F=\text{\rm tr} \left[\mathbf{M}\mathbf{M}^\top\right]\,,
#' } while (squared) Frobenius loss analogously is defined by
#' \deqn{
#' \code{LossFrobenius2}(\mathbf{M}, g(\boldsymbol{\Sigma}))= ||\mathbf{M}-g(\mathbf{\boldsymbol{\Sigma}})||^2_F\,,
#' } where \eqn{M} here denotes a suitable estimator and the function \eqn{g} is
#' \eqn{g(x)=x} when \code{type="covariance matrix"} is
#' chosen, otherwise \eqn{g(x)=1/x} for \code{type="precision matrix"}. In case 
#' \code{normalized=TRUE} (default) the above losses are normalized by the matrix
#' dimension \eqn{p}. Furthermore we present alternative loss measures focused directly
#' on the spectra of the matrices: \code{DistanceEuclideanEigenvalues2} and 
#' \code{LossEuclideanEigenvalues2} defined as follows
#'  \deqn{
#'  \code{DistanceEuclideanEigenvalues2}(\mathbf{M}_1, \mathbf{M}_2)=
#'  \sum\limits_{i=1}^{p} (\lambda_i(\mathbf{M}_1)-\lambda_i(\mathbf{M}_2))^2\,
#'  } and
#'  \deqn{
#'  \code{LossEuclideanEigenvalues2}(\mathbf{M}, g(\boldsymbol{\Sigma}) )= 
#'  \sum\limits_{i=1}^{p} (\lambda_i(\mathbf{M})-\lambda_i(g(\boldsymbol{\Sigma})))^2
#'  }
#' where \eqn{\lambda_i(\mathbf{A})} are the eigenvalues of a generic matrix
#'  \eqn{\mathbf{A}}. Similarly, \code{normalized=TRUE} (by default) normalizes 
#'  the losses and the function \eqn{g(x)} is either \eqn{x} or \eqn{1/x}, which 
#'  can be specified in  \code{type}. 
#'  The last loss function is specifically designed
#'  for the portfolios and is measuring the variance of the given portfolio relatively
#'   to the variance of the Global Minimum Variance (GMV) portfolio. For
#'  a generic \eqn{p}-dimensional vector of portfolio weights \eqn{\mathbf{w}}
#'   it is defined as
#'  \deqn{
#'  \code{LossRelativeOutOfSampleVariance}(\mathbf{w}, \boldsymbol{\Sigma})= 
#'  \frac{V_{\mathbf{w}} - V_{GMV}}{V_{GMV}}\,,
#'  }
#'  where \eqn{V_{\mathbf{w}}=\mathbf{w}^\top\boldsymbol{\Sigma}\mathbf{w}} and \eqn{V_{GMV}=
#'  (\mathbf{1}_p^\top\boldsymbol{\Sigma}^{-1}\mathbf{1}_p)^{-1}}
#'  are the variances of the given portfolio \eqn{\mathbf{w}} and of the GMV portfolio,
#'  respectively. The vector \eqn{\mathbf{1}_p} is the \eqn{p}-dimensional vector of ones.
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
#'   \code{NormFrobenius2(x \%*\% Sigma - diag(p) ) / p}.
#'   
#'   \item \code{type = "covariance matrix"} corresponds to the normalized
#'   Frobenius loss for the estimation of the precision matrix, i.e.
#'   \code{NormFrobenius2(x - Sigma ) / p}.
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
#' \code{NormFrobenius2} returns the squared Frobenius norm of a matrix.
#' \code{LossFrobenius2} returns the (normalized) Frobenius loss.
#'
#' \code{LossRelativeOutOfSampleVariance} returns a positive numeric value,
#' with attributes \code{"V_portfolio"} and \code{"V_GMV"}, which are respectively
#' the (out of sample) variances of the given \code{portfolioWeights} and of the
#' GMV portfolio.
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
#' 
#' LossFrobenius2(estimatedCov_sample, M, type = "covariance")
#' LossEuclideanEigenvalues2(estimatedCov_sample, M, type = "covariance")
#' 
#' 
#' # Losses for portfolios
#' Sigma = diag(1:5)
#' 
#' X <- MASS::mvrnorm(n = 100, mu = rep(0,5), Sigma = Sigma)
#' weights100 = GMV_Moore_Penrose(X)
#' 
#' X <- MASS::mvrnorm(n = 3, mu = rep(0,5), Sigma = Sigma)
#' weights3 = GMV_Moore_Penrose(X)
#' 
#' LossRelativeOutOfSampleVariance(weights100, Sigma)
#' LossRelativeOutOfSampleVariance(weights3, Sigma)
#' 
#' trueWeights = rowSums(solve(Sigma)) / sum(solve(Sigma))
#' 
#' NormFrobenius2(trueWeights - weights100, normalized = FALSE)
#' NormFrobenius2(trueWeights - weights3, normalized = FALSE)
#' 
#' @export
#' @rdname quadratic_losses
LossFrobenius2 <- function(x, Sigma, type, normalized = TRUE, ...) {
  UseMethod("LossFrobenius2")
}


#' @export
#' @rdname quadratic_losses
LossFrobenius2.matrix <- function(x,
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
      result = NormFrobenius2(x %*% Sigma - diag(p), normalized = normalized)
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
    x, Sigma, type = "precision matrix", normalized = TRUE, ...)
{
  type = match.arg(type)
  
  if (type != "precision matrix"){
    stop("Type is chosen to be ", type,
         " but x is of class 'EstimatedPrecisionMatrix'.")
  }
  result = LossFrobenius2(as.matrix(x), Sigma,
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
  result = LossFrobenius2(as.matrix(x), Sigma,
                          type = "covariance matrix", normalized = normalized)
  
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


#' @export
#' @rdname quadratic_losses
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

