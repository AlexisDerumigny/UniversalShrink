

#' Optimal portfolio weights based on a plug-in of the Moore-Penrose inverse
#' of the sample covariance matrix
#' 
#'
#' @param Y data matrix (rows are features, columns are observations).
#' TODO: transpose everything.
#' 
#' @export
GMV_Moore_Penrose <- function(Y, centeredCov = TRUE)
{
  iS_MP = Moore_Penrose(Y = Y, centeredCov = centeredCov)
  GMV_MP = GMV_PlugIn(iS_MP)
  
  return (GMV_MP)
}


