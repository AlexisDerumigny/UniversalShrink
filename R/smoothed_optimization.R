
optimization <- function(FUN, optimizationControls = NULL,
                         maximum = FALSE, verbose = 0, ...)
{
  method = optimizationControls$method
  
  if (method == "smoothed") {
    
    result = smoothed_optimization(
      FUN = FUN, grid = optimizationControls$grid, k = optimizationControls$k,
      verbose = verbose, maximum = maximum, ...)
    
  } else if (method == "optim with tan") {
    
    if (is.null(optimizationControls$lower)) {
      optimizationControls$lower <- 1/(10^6)
    }
    if (is.null(optimizationControls$upper)) {
      optimizationControls$upper <- pi/2 - optimizationControls$lower
    }
    if (is.null(optimizationControls$initialValue)) {
      optimizationControls$initialValue <- 1.5
    }
    
    FUN2 = function (u){ return (FUN(tan(u), ...))}
    
    result = stats::optim(
      par = optimizationControls$initialValue, fn = FUN2,
      lower = optimizationControls$lower, upper = optimizationControls$upper,
      method = "L-BFGS-B",
      control = list(fnscale = -maximum) # -1 for maximization
    )
    
    result$optimal_t = tan(result$par)
    
  } else if (method == "optimize") {
    
    if (is.null(optimizationControls$lower)) {
      optimizationControls$lower <- 1/(10^6)
    }
    if (is.null(optimizationControls$upper)) {
      optimizationControls$upper <- 50
    }
    
    result <- stats::optimize(f = FUN, maximum = maximum,
                              interval = c(optimizationControls$lower,
                                           optimizationControls$upper), ...)
    if (maximum) {
      result$optimal_t = result$maximum
    } else {
      result$optimal_t = result$minimum
    }
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
  if (verbose > 0) {
    pb_optimization = pbapply::startpb(max = length(grid))
  }
  
  for (i_t in 1:length(grid)) {
    vec_loss[i_t] = FUN(t = grid[i_t], ...)
    if (verbose > 0) {
      pbapply::setpb(pb_optimization, i_t)
    }
  }
  if (verbose > 0) {
    pbapply::closepb(pb_optimization)
  }
  
  # Smoothing
  fit.rollmedian = stats::runmed(vec_loss, k = k)
  
  # Finding the maximum or minimum as requested
  if (maximum) {
    optimal_t_unsmoothed = grid[which.max( vec_loss )]
    optimal_t_smoothed = grid[which.max( fit.rollmedian )]
  } else {
    optimal_t_unsmoothed = grid[which.min( vec_loss )]
    optimal_t_smoothed = grid[which.min( fit.rollmedian )]
  }
  
  result = list(
    optimal_t_unsmoothed = optimal_t_unsmoothed,
    optimal_t = optimal_t_smoothed,
    optimization_type = "smoothed",
    k = k,
    grid = grid,
    vec_loss = vec_loss,
    vec_loss_smoothed = fit.rollmedian,
    maximum = maximum
  )
  
  return (result)
}

# Compute a regular grid using the eigenvalues of S.
#
# The length of this grid is capped at max_length.
grid_optimization_default <- function(S, c_n, p, n, max_length, verbose = 0){
  if (is.null(max_length)) {
    max_length = 100
  }
  
  eigenvalues_S = eigen(S)$values
  if (p - 1 > n) {
    eigenvalues_S_ = eigenvalues_S[1:(n - 1)]
  } else {
    eigenvalues_S_ = eigenvalues_S
  }
  epsilon_steps = stats::median(abs(diff(eigenvalues_S_))) / (1 + sqrt(c_n))^2
  from = epsilon_steps / 2
  to = max(eigenvalues_S_)
  by = epsilon_steps
  grid_optim = seq(from = from, to = to, by = by)
  length_grid = length(grid_optim)
  
  if (verbose > 0) {
    cat("Computation of the default grid for optimization purpose:",
        "from =", from, " , to = ", to, " , by = ", by, 
        " , length = ", length_grid, "\n")
  }
  
  if (length_grid > max_length) {
    grid_optim = seq(from = from, to = to, length.out = max_length)
    
    if (verbose > 0) {
      cat("Grid length ", length_grid, " exceeds max_length ", max_length,
          ". Now changed to use a grid of length 'max_length'.\n")
    }
  }
  
  return (grid_optim)
}
