
optimization <- function(FUN, optimizationMethod,
                         grid = NULL, k = NULL,
                         initialValue = NULL, lower = NULL, upper = NULL,
                         maximum = FALSE, verbose, ...)
{
  if (optimizationMethod == "smoothed") {
    result = smoothed_optimization(
      FUN = FUN, grid = grid, k = k, verbose = verbose, maximum = maximum, ...)
    
  } else if (optimizationMethod == "optim with tan") {
    FUN2 = function (u){ return (FUN(tan(u), ...))}
    
    result = stats::optim(
      par = initialValue, fn = FUN2, lower = lower, upper = upper,
      method = "L-BFGS-B",
      control = list(fnscale = -maximum) # -1 for maximization
    )
    
    result$optimal_t = tan(result$par)
  }
  
  return (result)
}


smoothed_optimization <- function(FUN, grid, k, verbose, maximum, ...)
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
  
  # Smoothing
  fit.rollmedian = stats::runmed(vec_loss, k = k)
  
  # Finding the maximum or minimum as requested
  if (maximum) {
    optimal_t = grid[which.max( vec_loss )]
    optimal_t_smoothed = grid[which.max( fit.rollmedian )]
  } else {
    optimal_t = grid[which.min( vec_loss )]
    optimal_t_smoothed = grid[which.min( fit.rollmedian )]
  }
  
  result = list(
    optimal_t = optimal_t,
    optimal_t_smoothed = optimal_t_smoothed,
    optimization_type = "smoothed",
    k = k,
    grid = grid,
    vec_loss = vec_loss,
    maximum = maximum
  )
  
  return (result)
}


make_default_grid_optimization <- function(S, c_n){
  eigenvalues_S = eigen(S)$values
  epsilon_steps = median(abs(diff(eigenvalues_S))) / (1 + sqrt(c_n))^2
  grid_optim = seq(from = epsilon_steps / 2,
                   to = max(eigenvalues_S), 
                   by = epsilon_steps)
}
