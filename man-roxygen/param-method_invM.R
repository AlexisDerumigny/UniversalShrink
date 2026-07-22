
#' @param method_invM method for inverting the matrix \eqn{M}. Current
#' possible choices are \code{solve} (usual inverse) and \code{ginv}
#' (Moore-Penrose generalize inverse via \code{MASS::\link[MASS]{ginv}}).
#' The default is the faster and more precise method \code{recursive}.
