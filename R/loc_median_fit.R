


#' Estimate local median fit
#'
#' This function fits y based on x through a (weighted) median using
#' the \code{npoints/2} neighborhood.
#'
#' @param x,y the x and y coordinates of the points.
#' @param fraction,npoints the fraction / number of the points that are
#'   considered for each fit. `npoints` is the argument that is used in the
#'   end it is at least one. Default: `fraction = 0.1` and
#'   `npoints = length(x) * fraction`.
#' @param weighted a boolean that indicates if a weighted median is calculated.
#' @param ignore_zeros should the zeros be excluded from the fit
#'
#' @details
#' This function is low-level implementation detail and should usually not be
#' called by the user.
#'
#'
#' @seealso `locfit`: a package dedicated to local regression.
#'
#' @examples
#'   x <- runif(n = 1000, max = 4)
#'   y <- rpois(n = 1000, lambda = x * 10)
#'
#'   plot(x, y)
#'   fit <- loc_median_fit(x, y, fraction = 0.1)
#'   points(x, fit, col = "red")
#'
#'
#' @export
loc_median_fit <- function(x, y, fraction = 0.1, npoints = max(1, round(length(x) * fraction)),
                           weighted = TRUE, ignore_zeros = FALSE){
  # Make sure npoints is valid
  npoints <- max(1, npoints)
  npoints <- min(length(x), npoints)

  stopifnot(length(x) == length(y))
  if(length(x) == 0){
    return(numeric(0L))
  }
  stopifnot(npoints > 0 && npoints <= length(x))
  ordered_y <- y[order(x)]
  half_points <- floor(npoints/2)
  start <- half_points + 1
  end <- length(x) - half_points
  weights <- dnorm(seq(-3, 3, length.out = half_points * 2 + 1))

  if(end < start){
    if(weighted){
      wm <- matrixStats::weightedMedian(y, w =  dnorm(seq(-3, 3, length.out = length(x))))
      return(rep(wm, length(x)))
    }else{
      return(rep(median(y), length(x)))
    }
  }

  res <- rep(NA, length(x))
  idx <- start

  while(idx <= end){
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
  res[seq(length(x), end+1)] <- res[end]

  res[order(order(x))]
}


