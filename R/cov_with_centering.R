

#' Compute the covariance matrix with possible centering and the corresponding concentration ratio
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
#' \code{concentration_ratio} returns a numeric vector of length 1.
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
concentration_ratio <- function(n, p, centeredCov, verbose){
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

