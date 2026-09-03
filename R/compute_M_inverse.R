


#' Computing the inverse of the matrix M
#' 
#' @param all_tr they should be computed differently for covariance matrix
#' estimation and for precision matrix estimation.
#' It is a vector of size at least 2m
#' 
#' @noRd
compute_M_inverse <- function(m, all_tr0, all_tr, verbose = 0, mpfr, precBits){
  
  if (mpfr) {
    all_tr0 = Rmpfr::mpfr(all_tr0, precBits = precBits)
    all_tr  = Rmpfr::mpfr(all_tr , precBits = precBits)
    old_invM = methods::new("mpfrMatrix", all_tr0, Dim = c(1L, 1L))
    
    zero = Rmpfr::mpfr(0, precBits = precBits)
  } else {
    old_invM = matrix(all_tr0, ncol = 1, nrow = 1)
  }
  for (m_ in 1:m){
    if (mpfr) {
      invM = methods::new("mpfrMatrix",
                          rep(zero, (m_ + 1)^2), Dim_ = c(m_ + 1L, m_ + 1L))
      
      m_tilde = methods::new("mpfrMatrix",
                             all_tr[m_:(2 * m_ - 1)], Dim = c(m_, 1L))
    } else {
      
      invM = matrix(ncol = m_ + 1, nrow = m_ + 1)
      
      m_tilde = matrix(all_tr[m_:(2 * m_ - 1)], ncol = 1)
    }
    
    # Transposed version:
    m_tilde_t = t(m_tilde)
    
    vector_correction_term = old_invM %*% m_tilde
    
    xi_n_m = all_tr[2 * m_] - m_tilde_t %*% vector_correction_term
    xi_n_m = as.numeric(xi_n_m)
    
    corr_term = vector_correction_term %*% t(vector_correction_term) / xi_n_m
    
    invM[1:m_, 1:m_] = old_invM + corr_term
    invM[1:m_, m_+1] = - vector_correction_term / xi_n_m
    invM[m_+1, 1:m_] = - t(vector_correction_term) / xi_n_m
    invM[m_+1, m_+1] = 1 / xi_n_m
    
    old_invM = invM
    if (verbose > 0){
      cat("m =", m_, "M^{-1} =\n")
      print(invM)
    }
  }
  
  return (invM)
}

