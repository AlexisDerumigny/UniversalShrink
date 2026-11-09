

#' @param optimizationControls It describes the optimization method that it
#' used to choose \code{t} (if not given). The input \code{optimizationControls}
#' is either \code{NULL} (default optimization), or a list. This list can have
#' the following components:
#' 
#' - \code{method}: the name of the optimization method that is used.
#'   Currently, the following possibilities are implemented:
#'   * \code{"smoothed"}: this computes the loss on a grid, smooths it with a
#'   rolling median and minimizes this smoothed loss.
#'   
#'   * \code{"optimize"}: this uses the \code{stats::\link[stats]{optimize}}
#'   function.
#'   
#'   * \code{"optim with tan"}: uses \code{\link[stats]{optim}}
#'   with the change of variable \eqn{t = \tan(u)}, where
#'   \eqn{u = \arctan(t)}, to optimize over the interval
#'   \eqn{u \in [10^{-6}, \pi/2 - 10^{-6}]} by default.
#'   
#' - \code{k}: an integer giving the width of the window for median smoothing.
#' Used only for \code{method = "smoothed"}.
#'   
#' - \code{grid,max_length}: \code{grid} is an increasing (not necessarily
#' uniform) grid of positive reals at which to evaluate the loss. If \code{grid}
#' is \code{NULL} (the default), a regular grid is computed using the
#' eigenvalues of the sample covariance matrix.
#' To reduce computation time, the length of this grid is capped at
#' \code{max_length} (by default, 100). Use \code{max_length = Inf} to use the
#' ideal grid.
#' Used only for \code{method = "smoothed"}.
#' 
#' - \code{lower,upper}: lower and upper bound on the search interval.
#' Used only for \code{method = "optimize"} and \code{"optim with tan"}.
#' 
#' - \code{initialValue}: initial value of the search.
#' Used only for \code{"optim with tan"}.
#' 

