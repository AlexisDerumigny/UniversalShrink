
estimator_d0_thetaknown <- function(Ip, Sn, t, Theta){
  
  result = t * tr( solve(Sn + t * Ip) %*% Theta )
  return (result)
}

estimator_d1_thetaknown <- function(Ip, Sn, t, Theta){
  
  iS_ridge = solve(Sn + t * Ip)
  
  numerator = t * tr( iS_ridge %*% iS_ridge %*% Theta ) - t * estimator_d0_thetaknown(Ip, Sn, t, Theta)
  denominator = cn
  
  result = numerator / denominator
  
  return (result)
}

