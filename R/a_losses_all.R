
#' All quadratic losses
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

