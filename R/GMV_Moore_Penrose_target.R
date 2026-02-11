
#' First-order shrinkage of the Moore-Penrose portfolio towards a general target 
#' portfolio \eqn{\mathbf{b}}
#'
#' This function computes
#' \deqn{\hat{\alpha}^*\times \mathbf{w}_{MP} + (1 - \hat{\alpha}^*) \times \mathbf{b}}
#' where \eqn{\hat{\alpha}^*} is given by
#' \deqn{
#' \hat{\alpha}^*=
#' \frac{\mathbf{b}^\top\mathbf{S}_n\mathbf{b}-
#' \frac{\hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma} \right)}
#' {\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}}
#' {\mathbf{b}^\top\mathbf{S}_n\mathbf{b} - 
#' 2\frac{\hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma} \right)}
#' {\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
#' +\frac{\hat{d}_3\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
#' {\hat{d}_1^2\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}},
#' } where
#' \deqn{
#' \hat{d}_1\left(\mathbf{1}\mathbf{b}^\top\boldsymbol{\Sigma}\right)   = \frac{1}{\hat{v}(0)}
#' \left[\frac{1}{\hat{v}(0)}\left(1-\hat{d}_0(0,\mathbf{1}\mathbf{b}^\top) \right) 
#' - \hat{d}_1(\mathbf{1}\mathbf{b}^\top) \right],
#'}
#' with \eqn{\hat{v}(0)}, \eqn{\hat{d}_1\left( \mathbf{1}\mathbf{b}^\top\right)}, 
#' \eqn{\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
#' \eqn{d_3(\frac{\mathbf{1}\mathbf{1}^\top}{p})}, and
#'  \eqn{\hat{d}_0\left(0, \mathbf{1}\mathbf{b}^\top\right)} given in (8) and 
#'  in (S.41), (S.44), and (S.45)  from the supplement of Bodnar and Parolya (2026).
#'  The vector \eqn{\mathbf{w}_{MP}} are the optimal portfolio weights 
#'  estimated as the plug-in of the Moore-Penrose estimate of the precision matrix
#'   and \eqn{\mathbf{1}} is a vector of ones of size \eqn{p}.
#'
#'
# In the particular case of the shrinkage to the equally weighted portfolio,
# this function computes
# \deqn{\hat{\alpha}^*\times \mathbf{w}_{MP} + (1 - \hat{\alpha}^*) \times \mathbf{1}/p}
# where \eqn{\hat{\alpha}^*} is given by
# \deqn{
# \hat{\alpha}^*=
# \frac{\frac{1}{p}\mathbf{1}^\top\mathbf{S}_n\mathbf{1}-
# \frac{\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\boldsymbol{\Sigma} \right)}
# {\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}}
# {\frac{1}{p}\mathbf{1}^\top\mathbf{S}_n\mathbf{1} - 
# 2\frac{\hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\boldsymbol{\Sigma} \right)}
# {\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
# +\frac{\hat{d}_3\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}
# {\hat{d}_1^2\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}},
# } where
# \deqn{
# \hat{d}_1\left(\frac{\mathbf{1}\mathbf{1}^\top}{p}\boldsymbol{\Sigma}\right)   = \frac{1}{\hat{v}(0)}
# \left[\frac{1}{\hat{v}(0)}\left(1-\hat{d}_0(0,\frac{\mathbf{1}\mathbf{1}^\top}{p}) \right) 
# - \hat{d}_1(\frac{\mathbf{1}\mathbf{1}^\top}{p}) \right],
#}
# with \eqn{\hat{v}(0)}, \eqn{\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
# \eqn{\hat{d}_1\left( \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)}, 
# \eqn{d_3(\frac{\mathbf{1}\mathbf{1}^\top}{p})}, and
#  \eqn{\hat{d}_0\left(0, \frac{\mathbf{1}\mathbf{1}^\top}{p}\right)} given in (8) and 
#  in (S.41), (S.44), and (S.45)  from the supplement of Bodnar and Parolya (2026).
#  The vector \eqn{w_{MP}} are the optimal portfolio weights 
#  estimated as the plug-in of the Moore-Penrose estimate of the precision matrix
#   and \eqn{\mathbf{1}} is a vector of ones of size \eqn{p}.
#
#'
#'
#' @param X data matrix (rows are observations, columns are features).
#' 
#' @param b shrinkage target. By default, the equally-weighted portfolio is used
#' as a target.
#' 
#' @inheritParams cov_with_centering
#' 
#' @return a vector of size \eqn{p} of (estimated) optimal portfolio weights,
#' where \eqn{p} is the number of assets.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
#' 
#' @examples
#' set.seed(1)
#' n = 50
#' p = 2 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse (no shrinkage)
#' GMV_MP = GMV_Moore_Penrose(X)
#' 
#' Loss_GMV_MP = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_MP, Sigma = Sigma)
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse with shrinkage
#' # towards the equally weighted portfolio
#' GMV_MP_shrink_eq = GMV_Moore_Penrose_target(X)
#' 
#' Loss_GMV_MP_shrink_eq = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_MP_shrink_eq, Sigma = Sigma)
#' 
#' 
#' # Shrinkage helps to reduce the loss
#' stopifnot(Loss_GMV_MP_shrink_eq < Loss_GMV_MP)
#' 
#' 
#' # Compute GMV portfolio based on the Moore-Penrose inverse with shrinkage
#' # towards the true GMV portfolio
#' GMV_true = GMV_PlugIn(solve(Sigma))
#' 
#' GMV_MP_shrink_oracle = GMV_Moore_Penrose_target(X, b = GMV_true)
#' 
#' Loss_GMV_MP_shrink_oracle = LossRelativeOutOfSampleVariance(
#'   portfolioWeights = GMV_MP_shrink_oracle, Sigma = Sigma)
#' 
#' cat("GMV Moore-Penrose no shrinkage:" , Loss_GMV_MP, "\n")
#' cat("GMV Moore-Penrose target eq:"    , Loss_GMV_MP_shrink_eq, "\n")
#' cat("GMV Moore-Penrose_target oracle:", Loss_GMV_MP_shrink_oracle, "\n")
#' 
#' 
#' @export
GMV_Moore_Penrose_target <- function(X, centeredCov = TRUE, b = NULL,
                                     verbose = 0){
  if (is.null(b)){
    if (verbose > 0){
      cat("Default target: equally weighted portfolio\n")
    }
    result = GMV_Moore_Penrose_target_eq(
      X = X, centeredCov = centeredCov, verbose = verbose)
  } else {
    if (verbose > 0){
      cat("User-provided target portfolio\n")
    }
    result = GMV_Moore_Penrose_target_general(
      X = X, centeredCov = centeredCov, b = b, verbose = verbose)
  }
  
  return (result)
}

