

#' Compute the covariance matrix with possible centering and the corresponding 
#' concentration ratio
#' 
#' Having the number of dimensions \eqn{p}, the sample size \eqn{n} and the 
#' corresponding observation matrix \eqn{\mathbf{X}_n}, the user can choose 
#' between two versions of the sample covariance matrix: centered (mean is unknown) 
#' and noncentered one (mean is known to be zero). The noncentered sample covariance matrix 
#' (\code{centeredCov=FALSE}) is constructed simply by 
#' \deqn{\mathbf{S}_n=\frac{1}{n}\mathbf{X}_n^\top\mathbf{X}_n\,.
#' }  
#' In the case of unknown mean vector, one has to center the observation matrix
#' in the following way: 
#' let \eqn{\mathbf{J}_n = \mathbf{I}_n - \frac{1}{n}\mathbf{1}_n \mathbf{1}_n^\top}
#' with \eqn{\mathbf{1}} being the \eqn{n}-dimensional vector with all entries equal to one.
#' Then, \eqn{\mathbf{J}_n} is a projection matrix of rank \eqn{n-1}.
#' As such, all \eqn{n-1} nonzero eigenvalues of \eqn{\mathbf{J}_n} are equal to
#'  one, and its singular value decomposition is
#'\deqn{
#'  \mathbf{J}_n = \mathbf{H}_n \mathbf{H}_n^\top.
#'}
#'where \eqn{\mathbf{H}_n} is a \eqn{n \times (n-1)} orthogonal matrix, i.e.,
#'\eqn{\mathbf{H}_n^\top \mathbf{H}_n = \mathbf{I}_{n-1}}. Then, the centered sample
#'covariance matrix (\code{centeredCov=TRUE}) is defined in the following way
#'\deqn{
#'  \mathbf{S}_n
#'  = \frac{1}{n-1} \mathbf{X}^\top_n \mathbf{J}_n \mathbf{X}_n
#'  = \frac{1}{n-1} \tilde{\mathbf{X}}^\top_n \tilde{\mathbf{X}}_n,
#'  \quad \text{with} \quad
#'  \tilde{\mathbf{X}}_n = \mathbf{H}^\top_n\mathbf{X}_n.
#'  }
#'Similarly, the concentration ratio in the noncentered case is defined by 
#'\eqn{c_n=\frac{p}{n}}, while in the centered one it is given by 
#'\eqn{c_n=\frac{p}{n-1}}. Even though this difference looks somewhat subtle, it 
#'can have some significant influence on the performance of the estimators in 
#'finite sample regime.
#' 
#' @param X,n,p data matrix (rows are observations, columns are features).
#' \code{n} is the number of observations of the data matrix, and \code{p}
#' is the number of features
#' 
#' @param centeredCov Boolean variable. If \code{TRUE}, the covariance matrix is
#' computed using centering (i.e. in the general case where the mean of the
#' random vector may be non-zero).
#' If \code{FALSE} the covariance matrix is computed assuming that the mean of
#' the random vector of interest is \code{0} (i.e. does not need to be estimated).
#' 
#' @param verbose a number indicating whether to print intermediary values and
#' details about the progress of the computations. A value of \code{0} indicates
#' no printing at all, while higher values indicate increasingly more detailed
#' (more verbose) output. 
#' 
#' @returns \code{cov_with_centering} returns a \code{p * p} matrix.
#' \code{concentr_ratio} returns a numeric vector of length 1.
#' 
#' @export
#' 
cov_with_centering <- function(X, centeredCov){
  n = nrow(X)
  
  if (centeredCov){
    Jn <- diag(n) - matrix(1/n, nrow = n, ncol = n)
    
    # S <- Y %*% Jn %*% t(Y) / (n-1)
    S <- t(X) %*% Jn %*% X / (n-1)
  } else {
    
    # S <- Y %*% t(Y) / n
    S <- t(X) %*% X / n
  }
  
  return (S)
}


#' @rdname cov_with_centering
#' @export
concentr_ratio <- function(n, p, centeredCov, verbose){
  if (verbose > 0){
    cat("*  n = ", n, "\n")
    cat("*  p = ", p, "\n")
  }
  
  if (centeredCov){
    if (verbose > 0){
      cat("*  centered case\n")
    }
    
    c_n = p / (n-1)
    
  } else {
    if (verbose > 0){
      cat("*  non-centered case\n")
    }
    
    c_n = p / n
  }
  
  if (verbose > 0){
    cat("*  c_n = ", c_n, "\n\n")
  }
  
  return (c_n)
}

