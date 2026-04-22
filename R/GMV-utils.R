

#' Prepare and check the target vector in the GMV case
#' 
#' @param b target vector (vector of portfolio weights towards which shrinkage
#' is performed).
#' 
#' @param p dimension of the random vector of interest ( = number of assets
#' in the portfolio).
#' 
#' @returns a non-\code{NULL} vector b, checked.
#' 
#' @noRd
prepare_and_check_b <- function(b, p, tolerance = 0.001)
{
  if (is.null(b)){
    b = rep(1/p, length = p)
    
    return (b)
  }
  
  # Now we know that `b` is non-NULL, i.e. entered by the user.
  # so we try to convert it to a numeric vector.
  b = as.numeric(b)
  if (anyNA(b)){
    UniversalShrink_warning_condition_base(
      paste0(length(which(is.na(b))), " NAs in target portfolio 'b'."))
  }
  if (length(b) != p){
    stop("'b' should be a vector of length 'p'.")
  }
  if (abs(sum(b) - 1) > tolerance){
    stop("The weights (b) should sum up to 1.")
  }
  
  return (b)
}

