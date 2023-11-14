


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
#' @param sample_fraction fraction of the data to estimate local median on
#' @param seed control the random sampling of `sample_fraction`
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
#'   fit2 <- loc_median_fit(x, y, fraction = 0.1, sample_fraction = 0.75)
#'   points(x, fit, col = "red")
#'   points(x, fit, col = "blue")
#'
#'
#' @export
loc_median_fit <- function(x, y, fraction = 0.1, npoints = max(1, round(length(x) * fraction)),
                           weighted = TRUE, ignore_zeros = FALSE, sample_fraction = 1, seed=42){
  if (length(x) != length(y)) {
    stop("x and y must be of the same length.")
  }
  if (length(x) == 0) {
    return(numeric(0L))
  }
  if (!is.numeric(sample_fraction) || sample_fraction <= 0 || sample_fraction > 1) {
    stop("sample_fraction must be a numeric value between 0 and 1.")
  }
  if (is.null(npoints)) {
    npoints <- max(1, round(length(x) * fraction))
  }
  if (npoints <= 0 || npoints > length(x)) {
    stop("npoints must be a positive integer less than or equal to the length of x.")
  }

  if(sample_fraction != 1){
    set.seed(seed)
    subset_size <- round(sample_fraction * length(x))
    sample_indices <- sample(1:length(x), size = subset_size)
    x_sample <- x[sample_indices]
    y_sample <- y[sample_indices]
  } else {
    x_sample <- x
    y_sample <- y
  }
  
  # Make sure npoints is valid
  npoints <- max(1, npoints)
  npoints <- min(length(x_sample), npoints)

  ordered_indices <- order(x_sample)
  ordered_y <- y_sample[ordered_indices]
  half_points <- floor(npoints / 2)
  start <- half_points + 1
  end <- length(x_sample) - half_points
  weights <- dnorm(seq(-3, 3, length.out = half_points * 2 + 1))
  
  if(end < start){
    if(weighted){
      wm <- matrixStats::weightedMedian(y_sample, w =  dnorm(seq(-3, 3, length.out = length(x_sample))))
      return(rep(wm, length(x)))
    }else{
      return(rep(median(y_sample), length(x)))
    }
  }
  
  res <- rep(NA, length(x_sample))
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
  
  # Handle edge cases by extending the closest computed median
  res[seq_len(start - 1)] <- res[start]
  res[seq(end + 1, length(x_sample))] <- res[end]
  
  if (sample_fraction == 1) {
    # Reorder the results to match the original x
    res_ordered <- res
    res_ordered[ordered_indices] <- res
    return(res_ordered)
  } else {
    interp_func <- approxfun(x_sample[ordered_indices], res, method = "linear", yleft = res[start], yright = res[end])
    return(interp_func(x))
  }
}