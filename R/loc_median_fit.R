


#' Estimate local median fit
#'
#' This function fits y based on x through a (weighted) median using
#' the \code{npoints/2} neighborhood.
#'
#'
#'
#' @keywords internal
loc_median_fit <- function(x, y, fraction = 0.1, npoints = round(length(x) * fraction), weighted = TRUE, ignore_zeros = FALSE){
  stopifnot(length(x) == length(y))
  stopifnot(npoints > 0 && npoints <= length(x))
  ordered_y <- y[order(x)]
  half_points <- floor(npoints/2)
  start <- half_points + 1
  end <- length(x) - half_points

  res <- rep(NA, length(x))
  idx <- start
  weights <- dnorm(seq(-3, 3, length.out = half_points * 2 + 1))

  while(idx < end){
    selection <- ordered_y[seq(idx - half_points, idx + half_points)]
    if(ignore_zeros){
      if(weighted){
        used_weights <- weights[selection != 0]
        if(length(used_weights) == 0) used_weights <- 1
      }
      selection <- selection[selection != 0]
      if(length(selection) == 0) selection <- NA_real_
    }else{
      used_weights <- weights
    }
    if(!weighted){
      res[idx] <- median(selection)
    }else{
      res[idx] <- matrixStats::weightedMedian(selection, w = used_weights, interpolate = FALSE)
    }
    idx <- idx + 1
  }

  # Fill up NA's at the beginning and the end
  res[seq(1, start - 1)] <- res[start]
  res[seq(length(x), end)] <- res[end - 1]

  res[order(order(x))]
}


