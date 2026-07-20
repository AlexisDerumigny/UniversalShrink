
optimization <- function(FUN, smoothed, grid = NULL, k = NULL, verbose, ...)
{
  if (smoothed) {
    result = smoothed_optimization(
      FUN = FUN, grid = grid, k = k, verbose = verbose, ...)
  }
  
  return (result)
}


smoothed_optimization <- function(FUN, grid, k, verbose, ...)
{
  if (verbose > 0) {
    cat("Starting `smoothed_optimization` ... \n")
    cat("k = ", k, "\n")
  }
  
  n_t = length(grid)
  vec_loss = rep(NA_real_, length(grid))
  for (i_t in 1:length(grid)) {
    vec_loss[i_t] = FUN(t = grid[i_t], ...)
  }
  
  optimal_t = vec_t[which.min( vec_loss )]
  
  # Smoothing
  fit.rollmedian = stats::runmed(vec_loss, k = k, fill = NA)
  optimal_t_smoothed = vec_t[which.min( fit.rollmedian )]
  
  result = list(
    optimal_t = optimal_t,
    optimal_t_smoothed = optimal_t_smoothed,
    optimization_type = "smoothed",
    k = k,
    grid = grid,
    vec_loss = vec_loss
  )
  
  return (result)
}

