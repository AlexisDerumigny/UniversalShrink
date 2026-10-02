

# Trace function of a matrix
tr <- function(M){
  if (is.matrix(M)){
    result = sum(diag(M))
  } else {
    n = nrow(M)
    diagM = M[cbind(1:n, 1:n)]
    result = sum(diagM)
  }
  return (result)
}


format_ <- function(x, ...){
  if (is.numeric(x)){
    return (format(x, ...))
  }
  return (Rmpfr::formatMpfr(x, ...))
}


check_Rmpfr <- function (mpfr){
  if (!mpfr) {
    return (NULL) 
  }
  if (!requireNamespace("Rmpfr", quietly = TRUE)) {
    stop(UniversalShrink_error_condition_base(
      "Package \"Rmpfr\" must be installed to use the higher-precision option.",
      subclass = "MissingPackageError"
    ) )
  }
}


# Conversion to matrix and vectors  ============================================


#' Conversion of estimated matrices to matrix class and other generics
#' @name as.matrix.Estimator
#' 
#' @param x,object object to be converted, or for which ones want to extract the
#' coefficients
#' @param ... other arguments passed from methods, currently ignored.
#' 
#' @return `as.matrix()` returns the underlying (estimated) matrix.
#' `as.double()` and `as.numeric()` (which is an alias for it) 
#' returns the underlying (estimated) portfolio weights.
#' 
#' `coef()` returns a named numeric vector containing the estimator coefficients.
#' 
#' `get_t()` returns the value of the regularization parameter `t`
#' (used for ridge-type estimators), i.e. a numeric scalar.
#' 
NULL


#' @rdname as.matrix.Estimator
#' @export
as.matrix.PrecisionMatrix <- function(x, ...){
  if (length(x$matrix) == 0) {
    stop(
      UniversalShrink_error_condition_base(
        message = paste("Invalid x object:", 
                        "x must not have an empty precision matrix."),
        subclass = "InvalidArgumentError") )
  }
  return (x$matrix)
}


#' @rdname as.matrix.Estimator
#' @export
as.matrix.CovarianceMatrix <- function(x, ...){
  if (length(x$matrix) == 0) {
    stop(
      UniversalShrink_error_condition_base(
        message = paste("Invalid x object:", 
                        "x must not have an empty covariance matrix."),
        subclass = "InvalidArgumentError" ) )
  }
  return (x$matrix)
}


#' @rdname as.matrix.Estimator
#' @export
as.double.PortfolioWeights <- function(x, ...){
  return (x$portfolio_weights)
}


# Extracting coefficients  =====================================================


#' @rdname as.matrix.Estimator
#' @export
coef.EstimatedPrecisionMatrix <- function(object, ...)
{
  get_coefficient(object)
}

#' @rdname as.matrix.Estimator
#' @export
coef.EstimatedCovarianceMatrix <- function(object, ...)
{
  get_coefficient(object)
}

#' @rdname as.matrix.Estimator
#' @export
coef.EstimatedPortfolioWeights <- function(object, ...)
{
  get_coefficient(object)
}

get_coefficient <- function(object, ...)
{
  alpha <- if (!is.null(object$alpha_optimal)) {
    object$alpha_optimal
  } else {
    object$alpha
  }
  
  beta <- if (!is.null(object$beta_optimal)) {
    object$beta_optimal
  } else {
    object$beta
  }
  
  if (is.null(alpha) && is.null(beta)) {
    return (numeric(0))
  }
  
  if (is.null(alpha) && !is.null(beta)) {
    stop(UniversalShrink_error_condition_base(
      message = paste0(
        "A beta coefficient is stored in `object`, but no alpha ",
        "coefficient is available."
      ),
      subclass = c(
        "InvalidCoefficientRepresentationError",
        "InternalError"
      ),
      call = sys.call(-1),
      object = object
    ))
  }
  
  # Higher-order routines sometimes store alpha as a one-column matrix.
  alpha <- as.numeric(alpha)
  
  if (length(alpha) == 0L) {
    stop(UniversalShrink_error_condition_base(
      message = "The stored alpha coefficient must not be empty.",
      subclass = c(
        "InvalidCoefficientRepresentationError",
        "InternalError"
      ),
      call = sys.call(-1),
      object = object
    ))
  }
  
  if (is.null(beta)) {
    if (length(alpha) == 1L) {
      names(alpha) <- "alpha"
    } else {
      names(alpha) <- paste0("alpha_", seq_along(alpha) - 1L)
    }
    
    return(alpha)
  }
  
  beta <- as.numeric(beta)
  
  if (length(alpha) != 1L || length(beta) != 1L) {
    stop(UniversalShrink_error_condition_base(
      message = paste0(
        "When both alpha and beta are stored, each must be a numeric ",
        "value of length one. Here, length(alpha) = ", length(alpha),
        " and length(beta) = ", length(beta), "."
      ),
      subclass = c(
        "InvalidCoefficientRepresentationError",
        "InternalError"
      ),
      call = sys.call(-1),
      object = object,
      alpha = alpha,
      beta = beta
    ))
  }
  
  result = c(alpha = alpha, beta = beta)
  return (result)
}


#' @rdname as.matrix.Estimator
#' @export
get_t <- function(object, ...)
{
  t <- if (!is.null(object$t_optimal)) {
    object$t_optimal
  } else {
    object$t
  }
  
  t <- as.numeric(t)
  
  return (t)
}


# Constructors  ================================================================

#' Constructor for warning conditions of the package
#'
#' @noRd
UniversalShrink_warning_condition_base <- function(message, subclass = NULL,
                                                   call = sys.call(-1), ...) {
  # warningCondition() automatically adds 'warning' and 'condition' to the class
  return (
    warningCondition(
      message = message,
      class = c(subclass, "UniversalShrinkWarning"), # We add a base warning class
      call = call,
      ... # Allows for additional custom fields
    )
  )
}

#' Constructor for error conditions of the package
#'
#' @noRd
UniversalShrink_error_condition_base <- function(message, subclass = NULL,
                                                 call = sys.call(-1), ...) {
  # errorCondition() automatically adds 'error' and 'condition' to the class
  return (
    errorCondition(
      message = message,
      class = c(subclass, "UniversalShrinkError"), # Base error class
      call = call,
      ... # Allows for additional custom fields
    )
  )
}

