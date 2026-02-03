

#' Ridge with target set to the identity
#' 
#' This function computes
#' \deqn{\alpha \widehat{\Sigma^{-1}}^{ridge}_t + \beta I_p}
#' where \eqn{\widehat{\Sigma^{-1}}^{ridge}_t = (S + t I_p)^{-1}},
#' \eqn{S} is the sample covariance matrix,
#' \eqn{\alpha} and \eqn{\beta} are real-valued coefficients
#' and \eqn{I_p} is the identity matrix of size \eqn{p}.
#' 
#' The function `ridge_target_identity()` requires the specification of all
#' \eqn{t, \alpha, \beta}.
#' The function `ridge_target_identity_semioptimal()` only requires the
#' specification of \eqn{t} and compute (asymptotically) optimal choices of
#' \eqn{\alpha} and \eqn{\beta}.
#' Finally, the function `ridge_target_identity_optimal()` compute the (asymptotically)
#' optimal choice of \eqn{t, \alpha, \beta}.
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
#' Nestor Parolya & Taras Bodnar (2024).
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
#' precision_ridge_target_Cent = ridge_target(Y = t(X), centeredCov = TRUE)
#'     
#' precision_ridge_target_NoCent = ridge_target(Y = t(X), centeredCov = FALSE)
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
#'     ridge_target(Y = t(X), centeredCov = TRUE, Pi0 = solve(Sigma))
#' 
#' FrobeniusLoss2(precision_ridge_target_Cent_oracle, Sigma = Sigma)
#' 
#' 
#' @export
ridge_target <- function(Y, centeredCov = TRUE, Pi0 = NULL,
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
      
      none = ridge_target_identity(Y = Y, centeredCov = centeredCov,
                                   t = t, alpha = alpha, beta = beta,
                                   verbose = verbose),
      
      alpha_beta = ridge_target_identity_semioptimal(Y = Y, centeredCov = centeredCov,
                                                     t = t, verbose = verbose),
      
      all = ridge_target_identity_optimal(Y = Y, centeredCov = centeredCov,
                                          verbose = verbose, eps = eps, upp = upp,
                                          initialValue = initialValue)
    )
  } else {
    
    result = switch(
      optimizationType,
      
      none = ridge_target_general(Y = Y, centeredCov = centeredCov,
                                  t = t, alpha = alpha, beta = beta, Pi0 = Pi0,
                                  verbose = verbose),
      
      alpha_beta = ridge_target_general_semioptimal(Y = Y, centeredCov = centeredCov,
                                                    t = t, Pi0 = Pi0, verbose = verbose),
      
      all = ridge_target_general_optimal(Y = Y, centeredCov = centeredCov,
                                         Pi0 = Pi0,
                                         verbose = verbose, eps = eps, upp = upp,
                                         initialValue = initialValue)
    )
  }
  
  return (result)
}

