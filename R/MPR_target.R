


#' Moore-Penrose-Ridge with general target
#' 
#' This function computes
#' \deqn{\widehat{\Sigma^{-1}}^{ridge}_t
#'       = (S + t I_p)^{-1} - t * (S + t I_p)^{-2}},
#' where \eqn{S} is the sample covariance matrix and \eqn{t} is a given
#' parameter.
#' 
#' 
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @param t,alpha,beta,eps,upp,initialValue \code{t}, \code{alpha} and
#' \code{beta} are parameters of the estimation. In the optimized version,
#' the loss is optimized with respect to \eqn{u = arctan(t)} over the interval
#' \code{[eps, upp]}, and the optimizer starts at the \code{initialValue}.
#' 
#' @param Pi0 shrinkage target. By default it is the identity matrix.
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
#' Y = t(X)
#' 
#' # Estimation with default parameters (optimization) and identity target
#' precision_MPR_optimal = MPR_target(Y = Y)
#' 
#' cat("loss = ", FrobeniusLoss2(precision_MPR_optimal, Sigma = Sigma),
#'     ", t opt = ", precision_MPR_optimal$t_optimal, 
#'     ", alpha opt = ", precision_MPR_optimal$alpha_optimal,
#'     ", beta opt = ", precision_MPR_optimal$beta_optimal, "\n", sep = "")
#' 
#' # Estimation with default parameters (optimization) and oracle target
#' precision_MPR_optimal_oracle = MPR_target(
#'   Y = Y, Pi0 = solve(0.99 * Sigma + 0.01 * diag(nrow = p)) )
#'   
#' cat("loss = ", FrobeniusLoss2(precision_MPR_optimal_oracle, Sigma = Sigma),
#'     ", t = ", precision_MPR_optimal_oracle$t_optimal, 
#'     ", alpha opt = ", precision_MPR_optimal_oracle$alpha_optimal,
#'     ", beta opt = ", precision_MPR_optimal_oracle$beta_optimal, "\n", sep = "")
#' 
#' # Trying suboptimal alpha and beta
#' t_opt = precision_MPR_optimal$t_optimal
#' 
#' precision_MPR = MPR_target(Y = Y, t = t_opt, alpha = 1, beta = 0)
#' cat("loss = ", FrobeniusLoss2(precision_MPR, Sigma = Sigma))
#' 
#' # Trying suboptimal tm alpha and beta
#' precision_MPR = MPR_target(Y = Y, t = 1, alpha = 1, beta = 0)
#'    
#' cat("loss = ", FrobeniusLoss2(precision_MPR, Sigma = Sigma))
#' 
#' # Comparing with the non-shrinked version
#' precision_MPR_no_shrink = MPR_no_shrinkage(Y = Y, t = t_opt)
#'                                       
#' cat("loss = ", FrobeniusLoss2(precision_MPR_no_shrink, Sigma = Sigma))
#' 
#' 
#' @export
#' 
MPR_target <- function(Y, centeredCov = TRUE, Pi0 = NULL,
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
      
      none = MPR_target_identity(Y = Y, centeredCov = centeredCov,
                                 t = t, alpha = alpha, beta = beta,
                                 verbose = verbose),
      
      alpha_beta = MPR_target_identity_semioptimal(Y = Y, centeredCov = centeredCov,
                                                   t = t, verbose = verbose),
      
      all = MPR_target_identity_optimal(Y = Y, centeredCov = centeredCov,
                                        verbose = verbose, eps = eps, upp = upp,
                                        initialValue = initialValue)
    )
  } else {
    
    result = switch(
      optimizationType,
      
      none = MPR_target_general(Y = Y, centeredCov = centeredCov,
                                t = t, alpha = alpha, beta = beta, Pi0 = Pi0,
                                verbose = verbose),
      
      alpha_beta = MPR_target_general_semioptimal(Y = Y, centeredCov = centeredCov,
                                                  t = t, Pi0 = Pi0, verbose = verbose),
      
      all = MPR_target_general_optimal(Y = Y, centeredCov = centeredCov,
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

