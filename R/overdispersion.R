


#' Estimate the Overdispersion for a Vector of Counts
#'
#' @param y a numeric or integer vector or matrix with the counts for which
#'   the overdispersion is estimated
#' @param mean a numeric vector of either length 1 or `length(y)` or if
#'  `y` is a matrix, a matrix with the same dimensions. Contains
#'   the predicted value for that sample. If missing: `mean(y)` / `rowMeans(y)`
#' @param model_matrix a numeric matrix that specifies the experimental
#'   design. It can be produced using `stats::model.matrix()`.
#'   Default: `NULL`
#' @param do_cox_reid_adjustment the classical maximum likelihood estimator of the `overdisperion` is biased
#'   towards small values. McCarthy _et al._ (2012) showed that it is preferable to optimize the Cox-Reid
#'   adjusted profile likelihood.\cr
#'   `do_cox_reid_adjustment` can be either be `TRUE` or `FALSE` to indicate if the adjustment is
#'   added during the optimization of the `overdispersion` parameter. Default: `TRUE` if a
#'   model matrix is provided, otherwise `FALSE`
#' @inheritParams glm_gp
#'
#' @details
#' The function is optimized to be fast on many small counts. To
#' achieve this, the frequency table of the counts is calculated and
#' some parts are calculated based on this frequency table. If
#' there are probably many unique counts (heuristic: `max(y) > length(y)`),
#' the optimization is skipped.
#'
#' An earlier version of this package (< 1.1.1) used a separate
#' set of functions for the case of many small counts based on a paper
#' by Bandara et al. (2019). However, this
#' didn't bring a sufficient performance increase and meant an
#' additional maintenance burden.
#'
#' @return
#' The function returs a list with the following elements:
#' \describe{
#'   \item{`estimate`}{the numerical estimate of the overdispersion.}
#'   \item{`iterations`}{the number of iterations it took to calculate
#'   the result.}
#'   \item{`message`}{additional information about the fitting process.}
#' }
#'
#' @examples
#'  set.seed(1)
#'  # true overdispersion = 2.4
#'  y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
#'  # estimate = 1.7
#'  overdispersion_mle(y)
#'
#'
#'  # true overdispersion = 0
#'  y <- rpois(n = 10, lambda = 3)
#'  # estimate = 0
#'  overdispersion_mle(y)
#'  # with different mu, overdispersion estimate changes
#'  overdispersion_mle(y, mean = 15)
#'  # Cox-Reid adjustment changes the result
#'  overdispersion_mle(y, mean = 15, do_cox_reid_adjustment = FALSE)
#'
#'
#'  # Many very small counts, true overdispersion = 50
#'  y <- rnbinom(n = 1000, mu = 0.01, size = 1/50)
#'  summary(y)
#'  # estimate = 31
#'  overdispersion_mle(y, do_cox_reid_adjustment = TRUE)
#'
#'  # Function can also handle matrix input
#'  Y <- matrix(rnbinom(n = 10 * 3, mu = 4, size = 1/2.2), nrow = 10, ncol = 3)
#'  Y
#'  as.data.frame(overdispersion_mle(Y))
#'
#' @seealso [glm_gp()]
#' @export
overdispersion_mle <- function(y, mean,
                           model_matrix = NULL,
                           do_cox_reid_adjustment = ! is.null(model_matrix),
                           subsample = FALSE,
                           verbose = FALSE){

  # Validate n_subsampling
  stopifnot(length(subsample) == 1, subsample >= 0)
  if(isFALSE(subsample)){
    n_subsamples <- length(y)
  }else if(isTRUE(subsample)){
    n_subsamples <- 1000
  }else{
    n_subsamples <- subsample
  }

  if(! is.vector(y)){
    if(missing(mean)){
      mean <- array(DelayedMatrixStats::rowMeans2(y), dim = dim(y))
    }
    n_subsamples <- min(n_subsamples, ncol(y))
    if(n_subsamples < ncol(y)){
      if(verbose){ message("Subsample data to ", n_subsamples, " columns.") }
    }
    if(is.null(model_matrix)){
      model_matrix <- matrix(1, nrow = ncol(y), ncol = 1)
    }
    # This function calls overdispersion_mle() for each row, but is faster than a vapply()
    estimate_overdispersions_fast(y, mean, model_matrix, do_cox_reid_adjustment, n_subsamples)
  }else{
    overdispersion_mle_impl(as.numeric(y), mean, model_matrix, do_cox_reid_adjustment,
                            min(n_subsamples, length(y)), verbose = verbose)
  }

}



overdispersion_mle_impl <- function(y, mean,
                                 model_matrix,
                                 do_cox_reid_adjustment,
                                 n_subsamples,
                                 verbose = FALSE){

  stopifnot(is.numeric(y))
  if(missing(mean)){
    mean <- base::mean(y)
  }
  if(is.null(model_matrix)){
    model_matrix <- matrix(1, nrow = length(y), ncol = 1)
  }
  validate_model_matrix(model_matrix, matrix(y, nrow = 1))
  if(length(mean) == 1){
    mean <- rep(mean, length(y))
  }

  # Apply subsampling by randomly selecting elements of y and mean
  if(n_subsamples != length(y)){
    random_sel <- sort(sample(seq_along(y), size = round(n_subsamples), replace = FALSE))
    # It is important this is before subsetting y, because of lazy evaluation
    model_matrix <- model_matrix[random_sel, , drop=FALSE]
    y <- y[random_sel]
    mean <- mean[random_sel]
  }


  stopifnot(is.vector(y) && length(y) == length(mean))
  stopifnot(all(! is.na(y)))   # Cannot handle missing values
  stopifnot(all(y >= 0))
  stopifnot(all(! is.na(mean)))
  stopifnot(all(mean >= 0))
  stopifnot(all(is.finite(y)))
  stopifnot(all(is.finite(mean)))

  # Do conventional optimization
  conventional_overdispersion_mle(y, mean_vector = mean,
                                  model_matrix = model_matrix,
                                  do_cox_reid_adjustment = do_cox_reid_adjustment,
                                  verbose = verbose)
}



conventional_overdispersion_mle <- function(y, mean_vector,
                                       model_matrix = matrix(1, nrow = length(y), ncol = 1),
                                       do_cox_reid_adjustment = TRUE,
                                       verbose = FALSE){
  return_value = list(estimate = NA_real_, iterations = NA_real_, message = "")
  # Make a table of y
  if(length(y) == 0 || max(y) < length(y) * 10){
     tab <- make_table(y)
  }else{
     tab <- list(numeric(0), numeric(0))
  }

  if(all(y == 0)){
    return_value$estimate <- 0
    return_value$iterations <- 0
    return_value$message <- "All counts y are 0."
    return(return_value)
  }

  # Mu = 0 makes problems
  mean_vector[mean_vector == 0] <- 1e-6

  far_left_value <- conventional_score_function_fast(y, mu = mean_vector, log_theta = log(1e-8),
                                   model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
  if(far_left_value < 0){
    return_value$estimate <- 0
    return_value$iterations <- 0
    return_value$message <-  "Even for very small theta, no maximum identified"
    return(return_value)
  }

  mu <- mean(y)
  start_value <- (var(y) - mu) / mu^2
  if(is.na(start_value) || start_value <= 0){
    start_value <- 0.5
  }

  res <- nlminb(start = log(start_value),
         objective = function(log_theta){
           - conventional_loglikelihood_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
         }, gradient = function(log_theta){
           - conventional_score_function_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
         }, hessian = function(log_theta){
           res <- conventional_deriv_score_function_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
           matrix(- res, nrow = 1, ncol = 1)
         }, lower = log(1e-16), upper = log(1e16))

  if(res$convergence != 0){
    # Do the same thing again with numerical hessian as the
    # analytical hessian is sometimes less robust than the other
    # two functions
    res <- nlminb(start = log(start_value),
                  objective = function(log_theta){
                    - conventional_loglikelihood_fast(y, mu = mean_vector, log_theta = log_theta,
                                                      model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
                  }, gradient = function(log_theta){
                    - conventional_score_function_fast(y, mu = mean_vector, log_theta = log_theta,
                                                       model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment, tab[[1]], tab[[2]])
                  }, lower = log(1e-16), upper = log(1e16))
  }

  return_value$estimate <- exp(res$par)
  return_value$iterations <- res$iterations
  return_value$message <- res$message
  return_value
}





estimate_dispersions_by_moment <- function(Y, model_matrix, offset_matrix){
  xim <- 1/mean(DelayedMatrixStats::colMeans2(exp(offset_matrix)))
  bv <- DelayedMatrixStats::rowVars(Y)
  bm <- DelayedMatrixStats::rowMeans2(Y)
  (bv - xim * bm) / bm^2
}


estimate_dispersions_roughly <- function(Y, model_matrix, offset_matrix){
  # roughDisp <- DESeq2:::roughDispEstimate(y = Y / exp(offset_matrix),
  #                                         x = model_matrix)
  moments_disp <- estimate_dispersions_by_moment(Y, model_matrix, offset_matrix)
  # disp_rough <- pmin(roughDisp, moments_disp)
  disp_rough <- moments_disp
  ifelse(is.na(disp_rough) | disp_rough < 0, 0, disp_rough)
}





