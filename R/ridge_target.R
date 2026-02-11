

#' Ridge shrinkage estimator with general target
#' 
#' The optimal shrinkage estimator of the precision matrix with the ridge-type
#' inverse is given by
#' \deqn{
#' \widehat{\boldsymbol{\Pi}}_{R}
#' =
#' \hat{\alpha}_R^*(t^*)\,\mathbf{S}_n^{-}(t^*)
#' +
#' \hat{\beta}_R^*(t^*)\,\boldsymbol{\Pi}_0
#' }
#' where \eqn{\mathbf{S}_n^{-}(t^*) = (\mathbf{S}_n + t^*\mathbf{I}_p)^{-1}} and
#' \eqn{\mathbf{S}_n} is the sample covariance matrix. The coefficients
#' \eqn{\hat{\alpha}_R^*(t^*)} and \eqn{\hat{\beta}_R^*(t^*)} are given by
#' \deqn{
#' \hat{\alpha}_R^*(t^*)
#' =
#' \dfrac{
#' \hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' \hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#' }{
#' \left(
#' (t^*)^{-1}\hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' +
#' \hat{v}^{(1)}(t^*)\hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' (t^*)^{-1}\hat{d}_0^2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }
#' }
#' \deqn{
#' \hat{\beta}_R^*(t^*)
#' =
#' \dfrac{
#' \left(
#' (t^*)^{-1}\hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' +
#' \hat{v}^{(1)}(t^*)\hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \right)
#' \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#' -
#' (t^*)^{-1}
#' \hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }{
#' \left(
#' (t^*)^{-1}\hat{d}_0\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' +
#' \hat{v}^{(1)}(t^*)\hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' (t^*)^{-1}\hat{d}_0^2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }
#' }
#' The optimal penalty parameter is
#' \eqn{
#' t^* = \arg\max_{t>0}\, \hat{L}^2_{R;2}(t),
#' }
#' where
#' \deqn{
#' \hat{L}^2_{R;2}(t)
#' =
#' \frac{1}{\hat{q}_2\left(\dfrac{1}{p}\boldsymbol{\Pi}_0^2\right)}
#' \frac{
#' \left[
#' \hat{d}_0\left(t,\dfrac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{q}_2\left(\dfrac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' \hat{d}_0\left(t,\dfrac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' \hat{q}_1\left(\dfrac{1}{p}\boldsymbol{\Pi}_0\right)
#' \right]^2
#' }{
#' \left(
#' \hat{d}_0\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' +
#' t\,\hat{v}^{(1)}(t)\hat{d}_1\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' \hat{d}_0^2\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' .}
#' }
#' The quantities \eqn{\hat{v}^{(1)}(t)},
#' \eqn{\hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)},
#' \eqn{\hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)},
#' \eqn{\hat{d}_0\left(t,\frac{1}{p}\boldsymbol{\Sigma}\right)},
#' \eqn{\hat{d}_0\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)},
#' \eqn{\hat{d}_0\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)}
#' and \eqn{\hat{d}_1\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)} are defined
#' in the supplementary material of Bodnar and Parolya (2025). This procedure ensures
#'  that the loss 
#'\eqn{||\widehat{\boldsymbol{\Pi}}_{R}\boldsymbol{\Sigma}-\mathbf{I}_p||^2_F}
#' is asymptotically minimized with probability one.
#' 
#' Note that the function `ridge_target_identity()` requires the specification of all
#' \eqn{t, \alpha, \beta}.
#' The function `ridge_target_identity_semioptimal()` only requires the
#' specification of \eqn{t} and compute (asymptotically) optimal choices of
#' \eqn{\alpha} and \eqn{\beta}.
#' Finally, the function `ridge_target_identity_optimal()` compute the (asymptotically)
#' optimal choice of \eqn{t, \alpha, \beta}.
#' 
#' 
#' @param X data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param t,alpha,beta,eps,upp,initialValue \code{t}, \code{alpha} and
#' \code{beta} are parameters of the estimation. In the optimized version,
#' the loss is optimized with respect to \eqn{u = arctan(t)} over the interval
#' \code{[eps, upp]}, and the optimizer starts at the \code{initialValue}.
#' 
#' @param Pi0 shrinkage target. By default it is the identity matrix and
#' optimized computations are run in this case. Note that because of numerical
#' issues (with rounding), specifying explicitly \code{Pi0 = diag(p)} may give
#' different results for high values of \code{t}. This can arise as soon as
#' \code{t = 10^4}.
#' 
#' @inheritParams cov_with_centering
#' 
#' @returns the estimator of the precision matrix
#' (a `p` by `p` matrix).
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2026).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
#' 
#' 
#' @examples
#' 
#' n = 20
#' p = 5 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma=Sigma)
#' precision_ridge_target_Cent = ridge_target(X, centeredCov = TRUE)
#'     
#' precision_ridge_target_NoCent = ridge_target(X, centeredCov = FALSE)
#' 
#' FrobeniusLoss2(precision_ridge_target_Cent, Sigma = Sigma)
#' FrobeniusLoss2(precision_ridge_target_NoCent, Sigma = Sigma)
#' 
# # TODO: explore the convergence issues with good and bad target.
# # This works not as good as the MPR. The MPR takes a lot of time for the
# # (classical) optimization using the "L-BFGS-B" algorithm, but has finally
# # a smaller (true) loss. This is what we expected because we have a better
# # target.
#
#'
#' precision_ridge_target_Cent_oracle = 
#'   ridge_target(X, centeredCov = TRUE, Pi0 = solve(Sigma))
#' 
#' FrobeniusLoss2(precision_ridge_target_Cent_oracle, Sigma = Sigma)
#' 
#' 
#' @export
ridge_target <- function(X, centeredCov = TRUE, Pi0 = NULL,
                         t = NULL, alpha = NULL, beta = NULL,
                         verbose = 0,
                         eps = 1/(10^6), upp = pi/2 - eps, initialValue = 1.5)
{
  optimizationType = selectOptimizationType(t = t, alpha = alpha, beta = beta)
  
  if (verbose > 1){
    cat("optimizationType = '", optimizationType, "'.\n")
  }
  
  if (is.null(Pi0)){
    
    result = switch(
      optimizationType,
      
      none = ridge_target_identity(X = X, centeredCov = centeredCov,
                                   t = t, alpha = alpha, beta = beta,
                                   verbose = verbose),
      
      alpha_beta = ridge_target_identity_semioptimal(X = X, centeredCov = centeredCov,
                                                     t = t, verbose = verbose),
      
      all = ridge_target_identity_optimal(X = X, centeredCov = centeredCov,
                                          verbose = verbose, eps = eps, upp = upp,
                                          initialValue = initialValue)
    )
  } else {
    
    result = switch(
      optimizationType,
      
      none = ridge_target_general(X = X, centeredCov = centeredCov,
                                  t = t, alpha = alpha, beta = beta, Pi0 = Pi0,
                                  verbose = verbose),
      
      alpha_beta = ridge_target_general_semioptimal(X = X, centeredCov = centeredCov,
                                                    t = t, Pi0 = Pi0, verbose = verbose),
      
      all = ridge_target_general_optimal(X = X, centeredCov = centeredCov,
                                         Pi0 = Pi0,
                                         verbose = verbose, eps = eps, upp = upp,
                                         initialValue = initialValue)
    )
  }
  
  return (result)
}

