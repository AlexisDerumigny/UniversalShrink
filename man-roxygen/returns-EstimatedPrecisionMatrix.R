
#' @returns The estimated precision matrix, i.e. an estimator of the inverse
#' of the covariance matrix. It is an object of class
#' `c("EstimatedPrecisionMatrix", "PrecisionMatrix")`.
#' The underlying estimated precision matrix (of size \code{p * p}) can be
#' extracted using the method \code{\link[=as.matrix.Estimator]{as.matrix}}.
