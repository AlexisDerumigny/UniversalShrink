

#' Relative out-of-sample loss of a portfolio and Frobenius loss
#' 
#' The function \code{LossRelativeOutOfSampleVariance} returns the Variance of
#' the given portfolio relatively to the variance of the Global Minimum Variance
#' (GMV) portfolio.
#' For a generic \eqn{p}-dimensional vector of portfolio weights
#' \eqn{\mathbf{w}}, it is defined as
#' \deqn{
#'  \code{LossRelativeOutOfSampleVariance}(\mathbf{w}, \boldsymbol{\Sigma})= 
#'  \frac{V_{\mathbf{w}} - V_{GMV}}{V_{GMV}}\,,
#' }
#' where
#' \eqn{V_{\mathbf{w}} = \mathbf{w}^\top\boldsymbol{\Sigma}\mathbf{w}}
#' and
#' \eqn{V_{GMV} = (\mathbf{1}_p^\top \boldsymbol{\Sigma}^{-1} \mathbf{1}_p)^{-1}}
#' are the variances of the given portfolio \eqn{\mathbf{w}} and of the GMV
#' portfolio, respectively.
#' The vector \eqn{\mathbf{1}_p} is the \eqn{p}-dimensional vector of ones.
#' 
#' @param Sigma,SigmaInv the true covariance matrix and its inverse
#' 
#' @param x,portfolioWeights,otherPortfolioWeights vector of weights of a given
#' portfolio, of which we want to determine the loss (or an object of class
#' \code{EstimatedPortfolioWeights}).
#' 
#' @param normalized if \code{TRUE}, the Frobenius norm is divided by the number
#' of assets \code{p}.
#' 
#' @param ... Additional arguments passed to methods, currently ignored.
#' 
#' @returns \code{LossRelativeOutOfSampleVariance} returns a positive numeric
#' value, with attributes \code{"V_portfolio"} and \code{"V_GMV"}, which are
#' respectively the (out of sample) variances of the given
#' \code{portfolioWeights} and of the GMV portfolio.
#' 
#' \code{LossFrobenius2} returns a numeric value which is the Frobenius squared
#' norm between the two vectors of weights, see \code{\link{NormFrobenius2}}.
#' 
#' @examples
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
#' LossFrobenius2(weights100, trueWeights, normalized = FALSE)
#' LossFrobenius2(weights3, trueWeights, normalized = FALSE)
#' LossFrobenius2(weights100, weights3, normalized = FALSE)
#' 
#' @export
LossRelativeOutOfSampleVariance <- function(portfolioWeights, Sigma, SigmaInv = NULL){
  if (is.null(SigmaInv)){
    SigmaInv = solve(Sigma)
  }
  if (inherits(portfolioWeights, "EstimatedPortfolioWeights")){
    portfolioWeights = as.numeric(portfolioWeights)
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


#' @export
#' @rdname LossRelativeOutOfSampleVariance
LossFrobenius2.EstimatedPortfolioWeights <- function(
    x,
    otherPortfolioWeights,
    normalized = TRUE,
    ...)
{
  weights1 = as.numeric(x)
  weights2 = as.numeric(otherPortfolioWeights)
  
  result = NormFrobenius2(weights2 - weights1, normalized = normalized)
  return (result)
}


#' @export
#' @rdname LossRelativeOutOfSampleVariance
LossFrobenius2.numeric <- function(
    x,
    otherPortfolioWeights,
    normalized = TRUE,
    ...)
{
  weights1 = x
  weights2 = as.numeric(otherPortfolioWeights)
  
  result = NormFrobenius2(weights2 - weights1, normalized = normalized)
  return (result)
}
