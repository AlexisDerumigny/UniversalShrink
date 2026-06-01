


#' Computing the inverse of the matrix M
#' 
#' @param all_tr they should be computed differently for covariance matrix
#' estimation and for precision matrix estimation.
#' It is a vector of size at least 2m
#' 
#' @noRd
compute_M_inverse <- function(m, all_tr, p, verbose = 0){
  
  old_invM = matrix(1, ncol = 1, nrow = 1)
  for (m in 1:m){
    invM = matrix(ncol = m + 1, nrow = m + 1)
    
    m_tilde = matrix(all_tr[m:(2 * m - 1)] / p, ncol = 1)
    
    # Transposed version:
    m_tilde_t = t(m_tilde)
    
    vector_correction_term = old_invM %*% m_tilde
    
    xi_n_m = all_tr[2 * m] - m_tilde_t %*% vector_correction_term
    xi_n_m = as.numeric(xi_n_m)
    
    corr_term = vector_correction_term %*% t(vector_correction_term) / xi_n_m
    
    invM[1:m, 1:m] = old_invM + corr_term
    invM[1:m, m+1] = - vector_correction_term / xi_n_m
    invM[m+1, 1:m] = - t(vector_correction_term) / xi_n_m
    invM[m+1, m+1] = 1 / xi_n_m
    
    old_invM = invM
    if (verbose > 0){
      cat("m =", m, "M^{-1} =\n")
      print(invM)
    }
  }
  
  return (invM)
}

