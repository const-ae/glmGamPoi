



#' Estimate the Size Factors
#'
#' @param Y any matrix-like object (\code{base::matrix()}, \code{DelayedArray}, \code{HDF5Matrix},
#'   \code{Matrix::Matrix()}, etc.)
#'
#'
#'
#' Not exported
#' @keywords internal
estimate_size_factors <- function(Y){
  # Accept any matrix-like object
  log_geometric_means <- DelayedMatrixStats::rowMeans2(log(Y + 0.5))
  Y2 <- Y
  Y2[Y2 == 0] <- NA
  sf <- exp(DelayedMatrixStats::colMedians(subtract_vector_from_each_row(log(Y2), log_geometric_means), na.rm = TRUE))
  all_zero_column <- is.nan(sf)
  if(any(all_zero_column)){
    sf <- sf/exp(mean(log(sf), na.rm=TRUE))
    sf[all_zero_column] <- 0.001
  }else{
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}
