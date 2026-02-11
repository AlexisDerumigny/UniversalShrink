


#' Moore-Penrose-Ridge shrinkage with general target
#' 
#' Following Bodnar and Parolya (2006), the shrinkage estimator
#'  for the precision matrix using the Moore-Penrose-Ridge inverse of the sample 
#'  covariance matrix \eqn{\mathbf{S}_n} is computed as
#' \deqn{
#' \widehat{\boldsymbol{\Pi}}_{MPR}=\hat{\alpha}_{MPR}^*(t^*)\mathbf{S}_n^{\pm}(t^*)+
#' \hat{\beta}_{MP}^*(t^*)\boldsymbol{\Pi}_0\,,
#'} where \eqn{\mathbf{S}^+_n} denotes the Moore-Penrose inverse the sample
#' covariance matrix, \eqn{\boldsymbol{\Pi}_0} is the shrinkage target 
#' (\eqn{\boldsymbol{\Pi}_0=\mathbf{I}_p}, i.e., identity matrix by default) and
#'   \eqn{\hat{\alpha}_{MPR}^*(t^*)} and \eqn{\hat{\beta}_{MPR}^*(t^*)} are the optimal 
#'   shrinkage intensities given by 
#' \deqn{
#' \hat{\alpha}_{MPR}^*(t^*)
#' =
#' \dfrac{
#' -\hat{v}^{(1)}(t^*)
#' \hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' +
#' \hat{v}^{(1)}(t^*)
#' \hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#' }{
#' \hat{\grave{s}}_2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' [\hat{v}^{(1)}(t^*)]^2
#' \hat{d}_1^2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }
#' } and
#' \deqn{
#' \hat{\beta}_{MPR}^*(t^*)
#' =
#' \dfrac{
#' \hat{\grave{s}}_2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)
#' -
#' [\hat{v}^{(1)}(t^*)]^2
#' \hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{d}_1\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }{
#' \hat{\grave{s}}_2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' [\hat{v}^{(1)}(t^*)]^2
#' \hat{d}_1^2\left(t^*,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }.
#' }
#' To determine the optimal value of the tuning parameter we 
#' maximize \eqn{\hat{L}^2_{MPR;2}(t)} over \eqn{(0,\infty)} and obtain \eqn{t^*}. 
#' The function \eqn{\hat{L}^2_{MPR;2}(t)} is given by
#' \deqn{
#' \hat{L}^2_{MPR;2}(t)
#' =
#' \frac{1}{\hat{q}_2\left(\dfrac{1}{p}\boldsymbol{\Pi}_0^2\right)}
#' \frac{
#' [\hat{v}^{(1)}(t)]^2
#' \left[
#' \hat{d}_1\left(t,\dfrac{1}{p}\boldsymbol{\Sigma}\right)
#' \hat{q}_2\left(\dfrac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' \hat{d}_1\left(t,\dfrac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' \hat{q}_1\left(\dfrac{1}{p}\boldsymbol{\Pi}_0\right)
#' \right]^2
#' }{
#' \hat{\grave{s}}_2\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)
#' \hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)
#' -
#' [\hat{v}^{(1)}(t)]^2
#' \hat{d}_1^2\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)
#' }
#' }
#' where \eqn{\hat{v}^{(1)}(t)},
#' \eqn{\hat{q}_1\left(\frac{1}{p}\boldsymbol{\Pi}_0\right)},
#' \eqn{\hat{q}_2\left(\frac{1}{p}\boldsymbol{\Pi}_0^2\right)},
#' \eqn{\hat{\grave{s}}_2\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\right)},
#' \eqn{\hat{d}_1\left(t,\frac{1}{p}\boldsymbol{\Sigma}\right)}
#' and
#' \eqn{\hat{d}_1\left(t,\frac{1}{p}\boldsymbol{\Sigma}^2\boldsymbol{\Pi}_0\right)}
#' are defined in the supplementary material of Bodnar and Parolya (2026). 
#' This procedure ensures that the loss 
#'\eqn{||\widehat{\boldsymbol{\Pi}}_{MPR}\boldsymbol{\Sigma}-\mathbf{I}_p||^2_F}
#' is asymptotically minimized with probability one.
#' 
#' @param X data matrix (rows are observations, columns are features).
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
#' @returns the estimator of the precision matrix, of class
#' `EstimatedPrecisionMatrix`.
#' 
#' @references 
#' Nestor Parolya & Taras Bodnar (2024).
#' Reviving pseudo-inverses: Asymptotic properties of large dimensional
#' Moore-Penrose and Ridge-type inverses with applications.
#' \doi{10.48550/arXiv.2403.15792}
#' 
#' 
#' @examples
#' 
#' n = 10
#' p = 5 * n
#' mu = rep(0, p)
#' 
#' # Generate Sigma
#' X0 <- MASS::mvrnorm(n = 10*p, mu = mu, Sigma = diag(p))
#' H <- eigen(t(X0) %*% X0)$vectors
#' Sigma = H %*% diag(seq(1, 0.02, length.out = p)) %*% t(H)
#' 
#' # Generate example dataset
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' 
#' # Estimation with default parameters (optimization) and identity target
#' precision_MPR_optimal = MPR_target(X)
#' 
#' cat("loss = ", FrobeniusLoss2(precision_MPR_optimal, Sigma = Sigma),
#'     ", t opt = ", precision_MPR_optimal$t_optimal, 
#'     ", alpha opt = ", precision_MPR_optimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_optimal$beta_optimal, "\n", sep = "")
#' 
#' # Estimation with default parameters (optimization) and oracle target
#' oracle = solve(0.99 * Sigma + 0.01 * diag(nrow = p))
#' precision_MPR_optimal_oracle = MPR_target(X, Pi0 = oracle)
#'   
#' cat("loss = ", FrobeniusLoss2(precision_MPR_optimal_oracle, Sigma = Sigma),
#'     ", t = ", precision_MPR_optimal_oracle$t_optimal, 
#'     ", alpha opt = ", precision_MPR_optimal_oracle$alpha_optimal,
#'     ", beta opt = ", precision_MPR_optimal_oracle$beta_optimal, "\n", sep = "")
#' 
#' # Trying suboptimal alpha and beta
#' t_opt = precision_MPR_optimal$t_optimal
#' 
#' precision_MPR = MPR_target(X, t = t_opt, alpha = 1, beta = 0)
#' cat("loss = ", FrobeniusLoss2(precision_MPR, Sigma = Sigma))
#' 
#' # Trying suboptimal t, alpha and beta
#' precision_MPR = MPR_target(X, t = 1, alpha = 1, beta = 0)
#'    
#' cat("loss = ", FrobeniusLoss2(precision_MPR, Sigma = Sigma))
#' 
#' # Comparing with the non-shrinked version
#' precision_MPR_no_shrink = MPR_no_shrinkage(X, t = t_opt)
#'                                       
#' cat("loss = ", FrobeniusLoss2(precision_MPR_no_shrink, Sigma = Sigma))
#' 
#' 
#' @export
#' 
MPR_target <- function(X, centeredCov = TRUE, Pi0 = NULL,
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
      
      none = MPR_target_identity(X = X, centeredCov = centeredCov,
                                 t = t, alpha = alpha, beta = beta,
                                 verbose = verbose),
      
      alpha_beta = MPR_target_identity_semioptimal(X = X, centeredCov = centeredCov,
                                                   t = t, verbose = verbose),
      
      all = MPR_target_identity_optimal(X = X, centeredCov = centeredCov,
                                        verbose = verbose, eps = eps, upp = upp,
                                        initialValue = initialValue)
    )
  } else {
    
    result = switch(
      optimizationType,
      
      none = MPR_target_general(X = X, centeredCov = centeredCov,
                                t = t, alpha = alpha, beta = beta, Pi0 = Pi0,
                                verbose = verbose),
      
      alpha_beta = MPR_target_general_semioptimal(X = X, centeredCov = centeredCov,
                                                  t = t, Pi0 = Pi0, verbose = verbose),
      
      all = MPR_target_general_optimal(X = X, centeredCov = centeredCov,
                                       Pi0 = Pi0,
                                       verbose = verbose, eps = eps, upp = upp,
                                       initialValue = initialValue)
    )
  }
  
  return (result)
}



selectOptimizationType <- function (t, alpha, beta)
{
  if (is.null(t)){
    # if t is not given, we need full optimization anyway
    
    if (!is.null(alpha) || !is.null(beta)){
      warning(UniversalShrink_warning_condition_base(
        message = paste0(
          "The input t is NULL: ", t,
          " , but (at least) one of alpha and/or beta is provided:",
          "alpha = ", alpha, ", beta = ", "beta. \n",
          "Therefore optimization in t is done, and alpha and/or beta are ignored.",
          "If you want t alpha and beta to be used, please also specify t."
        ), subclass = "MissingParametersWarning") )
      
      optimizationType = "all"
    } else {
      # Nothing is specified so we need full optimization. This is the default case.
      optimizationType = "all"
    }
  } else {
    # t is specified. So we now look at alpha and beta
    
    if (is.null(alpha) && is.null(beta)) {
      optimizationType = "alpha_beta"
    } else if (is.null(beta)) {
      # alpha is not null but beta is null
      
      warning(UniversalShrink_warning_condition_base(
        message = paste0(
          "The input alpha is non-NULL: ",
          "alpha = ", alpha, ", but beta is NULL. \n",
          "Therefore alpha is ignored. If you want alpha to be used, please also ",
          "specify beta."
        ), subclass = "MissingParametersWarning") )
      
      optimizationType = "alpha_beta"
    } else if (is.null(alpha)) {
      # beta is not null but alpha is null
      
      warning(UniversalShrink_warning_condition_base(
        message = paste0(
          "The input beta is non-NULL: ",
          "beta = ", beta, ", but alpha is NULL. \n",
          "Therefore beta is ignored. If you want beta to be used, please also ",
          "specify alpha."
        ), subclass = "MissingParametersWarning") )
      
      optimizationType = "alpha_beta"
    } else {
      # all of t, alpha and beta are not NULL, so no optimization is needed.
      
      optimizationType = "none"
    }
  }
  
  return (optimizationType)
}

