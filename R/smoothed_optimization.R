

smoothed_optimization <- function(FUN, grid, k, verbose, ...)
{
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
    optimal_t_smoothed = optimal_t_smoothed
  )
  
  return (result)
}

