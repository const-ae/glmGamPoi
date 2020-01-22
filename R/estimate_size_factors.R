



#' Estimate the Size Factors
#'
#'
#'
#' Not exported
#' @keywords internal
estimate_size_factors <- function(Y){
  stopifnot(is.matrix(Y))
  logGeoMeans <- rowMeans(log(Y + 0.5))
  logGeoMeans_mat <- matrix(logGeoMeans, nrow = length(logGeoMeans), ncol=ncol(Y))
  Y2 <- Y
  Y2[Y2 == 0] <- NA
  sf <- exp(matrixStats::colMedians(log(Y2) - logGeoMeans_mat, na.rm=TRUE))
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
