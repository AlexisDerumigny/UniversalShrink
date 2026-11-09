
#' @param method_invM Character string specifying the method used for the
#' inversion of the matrix \eqn{M}. Possible choices are:
#' \itemize{
#'   \item \code{"recursive"}: compute the inverse using the recursive
#'     formula. This is the default and is generally faster and more
#'     numerically stable for larger values of \code{m}.
#'     
#'   \item \code{"solve"}: construct \eqn{M} and compute its usual inverse
#'     using \code{\link[base]{solve}}.
#'     
#'   \item \code{"ginv"}: construct \eqn{M} and compute its Moore-Penrose
#'     generalized inverse using \code{MASS::\link[MASS]{ginv}}.
#' }
