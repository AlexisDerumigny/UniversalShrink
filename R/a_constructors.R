# Internal constructors for estimator objects
#
# These constructors define the common representation of estimated precision
# matrices, estimated covariance matrices, and estimated portfolio weights.
#
# They are intentionally not exported.

stop_invalid_estimator_object <- function(message, estimate, call = sys.call(-1))
{
  stop(UniversalShrink_error_condition_base(
    message = message,
    subclass = c("InvalidEstimatorObjectError", "InternalError"),
    call = call,
    estimate = estimate
  ) )
}


check_estimated_matrix <- function(estimate, estimator_type, call = sys.call(-1))
{
  if (!is.matrix(estimate)) {
    stop_invalid_estimator_object(
      message = paste0("The estimated ", estimator_type, " must be a matrix."),
      estimate = estimate,
      call = call)
  }
  
  if (!is.numeric(estimate)) {
    stop_invalid_estimator_object(
      message = paste0(
        "The estimated ", estimator_type, " must be a numeric matrix."),
      estimate = estimate,
      call = call
    )
  }
  
  if (length(estimate) == 0L) {
    stop_invalid_estimator_object(
      message = paste0("The estimated ", estimator_type, " must not be empty."),
      estimate = estimate,
      call = call
    )
  }
  
  if (nrow(estimate) != ncol(estimate)) {
    stop_invalid_estimator_object(
      message = paste0(
        "The estimated ", estimator_type, " must be a square matrix. ",
        "Here, dim(estimate) = c(", paste(dim(estimate), collapse = ", "), ")."
      ),
      estimate = estimate,
      call = call
    )
  }
  
  invisible(NULL)
}


check_estimator_dimension <- function(p, expected_p, estimate, 
                                      call = sys.call(-1) )
{
  valid_p <- (is.numeric(p) && length(p) == 1L && !is.na(p) && is.finite(p) &&
      p == expected_p)
  
  if (!valid_p) {
    stop(UniversalShrink_error_condition_base(
      message = paste0(
        "`p` must equal the dimension of the estimate. ",
        "Here, p = ",
        paste(utils::capture.output(dput(p)), collapse = "\n"),
        " and the estimate has dimension ",
        expected_p,
        "."
      ),
      subclass = c("InvalidEstimatorObjectError", "InternalError"),
      call = call,
      estimate = estimate,
      p = p,
      expected_p = expected_p) )
  }
  
  invisible(NULL)
}


#' Construct an estimated precision matrix
#'
#' @param estimate Numeric square matrix containing the estimated precision
#'   matrix.
#' @param n Sample size, or `NA_integer_` when no sample size is associated
#'   with the estimate.
#' @param p Matrix dimension. It must agree with the dimension of `estimate`.
#' @param centeredCov A logical indicating whether centered covariance
#'   estimation was used, or `NA` when this is not applicable or unknown.
#' @param method Character description of the estimator.
#' @param call Matched call associated with the estimator.
#' @param ... Additional estimator-specific fields.
#'
#' @return An object of class `EstimatedPrecisionMatrix`.
#'
#' @noRd
new_estimated_precision_matrix <- function(
  estimate, n, p, centeredCov, method, call, ...)
{
  constructor_call <- sys.call()
  
  check_estimated_matrix(estimate = estimate,
                         estimator_type = "precision matrix",
                         call = constructor_call)
  
  check_estimator_dimension(p = p, expected_p = nrow(estimate), 
                            estimate = estimate, call = constructor_call)
  
  dots <- list(...)
  
  result <- c(list(estimated_precision_matrix = estimate,
                   n = n,
                   p = p,
                   centeredCov = centeredCov,
                   method = method,
                   call = call), 
              dots)
  
  class(result) <- "EstimatedPrecisionMatrix"
  
  return (result)
}


#' Construct an estimated covariance matrix
#'
#' @param estimate Numeric square matrix containing the estimated covariance
#'   matrix.
#' @inheritParams new_estimated_precision_matrix
#'
#' @return An object of class `EstimatedCovarianceMatrix`.
#'
#' @noRd
new_estimated_covariance_matrix <- function(
  estimate, n, p, centeredCov, method, call, ...)
{
  constructor_call <- sys.call()
  
  check_estimated_matrix(estimate = estimate,
                         estimator_type = "covariance matrix",
                         call = constructor_call)
  
  check_estimator_dimension(p = p, expected_p = nrow(estimate),
                            estimate = estimate, call = constructor_call)
  
  dots <- list(...)
  
  result <- c(list(estimated_covariance_matrix = estimate,
                   n = n,
                   p = p,
                   centeredCov = centeredCov,
                   method = method,
                   call = call),
              dots)
  
  class(result) <- "EstimatedCovarianceMatrix"
  
  return (result)
}


#' Construct estimated portfolio weights
#'
#' @param estimate Numeric vector containing the estimated portfolio weights.
#' @inheritParams new_estimated_precision_matrix
#'
#' @return An object of class `EstimatedPortfolioWeights`.
#'
#' @noRd
new_estimated_portfolio_weights <- function(
  estimate, n, p, centeredCov, method, call, ...)
{
  constructor_call <- sys.call()
  
  if (!is.numeric(estimate) || is.matrix(estimate) || is.array(estimate)) {
    stop_invalid_estimator_object(
      message = "The estimated portfolio weights must be a numeric vector.",
      estimate = estimate,
      call = constructor_call
    )
  }
  
  if (length(estimate) == 0L) {
    stop_invalid_estimator_object(
      message = "The estimated portfolio weights must not be empty.",
      estimate = estimate,
      call = constructor_call
    )
  }
  
  check_estimator_dimension(p = p, expected_p = length(estimate),
                            estimate = estimate, call = constructor_call)
  
  dots <- list(...)
  
  result <- c(list(estimated_portfolio_weights = estimate,
                   n = n,
                   p = p,
                   centeredCov = centeredCov,
                   method = method,
                   call = call),
              dots)
  
  class(result) <- "EstimatedPortfolioWeights"
  
  return (result)
}

