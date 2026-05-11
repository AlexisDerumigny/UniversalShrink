
#' All quadratic losses
#' 
#' Returns and print all (currently implemented) quadratic losses.
#' 
#' @param x an object, of class \code{EstimatedPrecisionMatrix},
#' \code{EstimatedCovarianceMatrix} or \code{AllLosses}.
#' 
#' @param Sigma,SigmaInv true covariance matrix and its inverse
#' (true precision matrix). If \code{SigmaInv} is needed and missing, it is
#' computed by numerical inversion of the provided \code{Sigma}.
#' @param ... other arguments, ignored except for the \code{print} method, for
#' which they are passed to \code{\link[base]{print.data.frame}}.
#' 
#' @returns \code{Losses} returns an object of class \code{AllLosses}.
#' \code{print} returns \code{NULL} and is only called for its side-effects.
#' 
#' \code{as.matrix}, \code{as.data.frame} and
#' \code{as.numeric} (or \code{as.double}) respectively
#' return a matrix, a data.frame or a numeric vector containing the losses.
#' They are all named.
#' 
#' @seealso The quadratic losses \code{\link{LossFrobenius2}},
#' \code{\link{LossInverseFrobenius2}} and \code{\link{LossEuclideanEigenvalues2}}.
#' 
#' 
#' @examples
#' 
#' n = 100
#' p = 200
#' Sigma = diag(seq(1, 0.02, length.out = p))
#' mu = rep(0, p)
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' 
#' precision_MP = Moore_Penrose(X)
#' Losses(precision_MP, Sigma = Sigma)
#' 
#' precision_MP_shrink = Moore_Penrose_shrinkage(X)
#' Losses(precision_MP_shrink, Sigma = Sigma)
#' 
#' estimatedCov_shrink = cov_analytical_NL_shrinkage(X)
#' Losses(estimatedCov_shrink, Sigma = Sigma)
#' 
#' L = Losses(estimatedCov_shrink, Sigma = Sigma)
#' as.matrix(L)
#' as.data.frame(L)
#' as.numeric(L)
#' as.double(L)
#' 
#' @export
Losses <- function(x, Sigma, ...) {
  UseMethod("Losses")
}

#' @rdname Losses
#' @export
Losses.EstimatedPrecisionMatrix <- function(x, Sigma, SigmaInv = NULL, ...)
{
  if (is.null(SigmaInv)) {
    SigmaInv = solve(Sigma)
  }
  Frob2 = c(LossFrobenius2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                           normalized = TRUE),
            LossFrobenius2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                           normalized = FALSE))
  
  InvFrob2 = c(
    LossInverseFrobenius2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                          normalized = TRUE),
    LossInverseFrobenius2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                          normalized = FALSE))
  
  EuclEig = c(
    LossEuclideanEigenvalues2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                              normalized = TRUE),
    LossEuclideanEigenvalues2(x = x, Sigma = Sigma, SigmaInv = SigmaInv,
                              normalized = FALSE))
  
  allLosses = rbind(Frobenius2 = Frob2,
                    InverseFrobenius2 = InvFrob2,
                    EuclideanEigenvalues = EuclEig)
  
  allLosses = as.data.frame(allLosses)
  colnames(allLosses) <- c("Normalized", "Unnormalized")
  
  result = list(
    allLosses = allLosses,
    estimated = x,
    problemType = "Estimation of the precision matrix"
  )
  class(result) <- "AllLosses"
  
  return (result)
}


#' @rdname Losses
#' @export
Losses.EstimatedCovarianceMatrix <- function(x, Sigma, ...)
{
  Frob2 = c(LossFrobenius2(x = x, Sigma = Sigma, normalized = TRUE),
            LossFrobenius2(x = x, Sigma = Sigma, normalized = FALSE))
  
  EuclEig = c(
    LossEuclideanEigenvalues2(x = x, Sigma = Sigma, normalized = TRUE),
    LossEuclideanEigenvalues2(x = x, Sigma = Sigma, normalized = FALSE))
  
  allLosses = rbind(Frobenius2 = Frob2,
                    EuclideanEigenvalues = EuclEig)
  
  allLosses = as.data.frame(allLosses)
  colnames(allLosses) <- c("Normalized", "Unnormalized")
  
  result = list(
    allLosses = allLosses,
    estimated = x,
    problemType = "Estimation of the covariance matrix"
  )
  class(result) <- "AllLosses"
  
  return (result)
}


#' @rdname Losses
#' @export
Losses.EstimatedPortfolioWeights <- function(x, Sigma, SigmaInv = NULL, ...)
{
  if (is.null(SigmaInv)){
    SigmaInv = solve(Sigma)
  }
  
  LossOutOfSampleVariance = c( 
    LossOutOfSampleVariance(portfolioWeights = x, Sigma = Sigma, 
                            SigmaInv = SigmaInv, normalized = TRUE),
    LossOutOfSampleVariance(portfolioWeights = x, Sigma = Sigma,
                            SigmaInv = SigmaInv, normalized = FALSE))
  
  GMV_portfolio = GMV_PlugIn(SigmaInv)
  
  
  Frob2 = c(LossFrobenius2(x = x, GMV_portfolio, normalized = TRUE),
            LossFrobenius2(x = x, GMV_portfolio, normalized = FALSE))
  
  allLosses = rbind(LossOutOfSampleVariance = LossOutOfSampleVariance,
                    Frobenius2 = Frob2)
  
  allLosses = as.data.frame(allLosses)
  colnames(allLosses) <- c("Normalized", "Unnormalized")
  
  result = list(
    allLosses = allLosses,
    estimated = x,
    problemType = "Estimation of portfolio weights"
  )
  class(result) <- "AllLosses"
  
  return (result)
}


#' @rdname Losses
#' @export
print.AllLosses <- function(x, ...){
  cat(x$problemType)
  cat(", method =", x$estimated$method)
  if (x$estimated$centeredCov){
    cat(" (centered)")
  } else {
    cat(" (non centered)")
  }
  cat("\n")
  cat("n =", x$estimated$n)
  cat(", p =", x$estimated$p)
  c_n = concentr_ratio(n = x$estimated$n,
                       p = x$estimated$p,
                       centeredCov = x$estimated$centeredCov, verbose = 0)
  cat(", c_n =", c_n)
  cat("\n\n")
  cat("Losses:\n")
  print(x$allLosses, ...)
  return (invisible(NULL))
}


#' @rdname Losses
#' @export
as.double.AllLosses <- function(x, ...){
  result = c(x$allLosses[,1], x$allLosses[,2])
  names(result) <- paste0(rep(rownames(x$allLosses), 2), "_",
                          rep(colnames(x$allLosses), each = nrow(x$allLosses)))
  return (result)
}


#' @rdname Losses
#' @export
as.data.frame.AllLosses <- function(x, ...){
  
  return (x$allLosses)
}


#' @rdname Losses
#' @export
as.matrix.AllLosses <- function(x, ...){
  result = as.matrix(x$allLosses)
  return (result)
}

