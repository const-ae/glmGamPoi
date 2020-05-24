


#' Estimate the Overdispersion for a Vector of Counts
#'
#' @param y a numeric or integer vector with the counts for which
#'   the overdispersion is estimated
#' @param mean a numeric vector of either length 1 or `length(y)`
#'   with the predicted value for that sample. Default: `mean(y)`.
#' @param model_matrix a numeric matrix that specifies the experimental
#'   design. It can be produced using `stats::model.matrix()`.
#'   Default: `matrix(1, nrow = length(y), ncol = 1)`, which is the
#'   model matrix for a 'just-intercept-model'.
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
#'  gampoi_overdispersion_mle(y)
#'
#'
#'  # true overdispersion = 0
#'  y <- rpois(n = 10, lambda = 3)
#'  # estimate = 0
#'  gampoi_overdispersion_mle(y)
#'  # with different mu, overdispersion estimate changes
#'  gampoi_overdispersion_mle(y, mean = 15)
#'  # Cox-Reid adjustment changes the result
#'  gampoi_overdispersion_mle(y, mean = 15, do_cox_reid_adjustment = FALSE)
#'
#'
#'  # Many very small counts, true overdispersion = 50
#'  y <- rnbinom(n = 1000, mu = 0.01, size = 1/50)
#'  summary(y)
#'  # estimate = 31
#'  gampoi_overdispersion_mle(y)
#'
#' @seealso [glm_gp()]
#' @export
gampoi_overdispersion_mle <- function(y, mean = base::mean(y),
                           model_matrix = matrix(1, nrow = length(y), ncol = 1),
                           do_cox_reid_adjustment = TRUE,
                           subsample = FALSE,
                           verbose = FALSE){
  y <- as.numeric(y)
  validate_model_matrix(model_matrix, matrix(y, nrow = 1))
  if(length(mean) == 1){
    mean <- rep(mean, length(y))
  }
  # Validate n_subsampling
  stopifnot(length(subsample) == 1, subsample >= 0)
  if(isFALSE(subsample)){
    n_subsamples <- length(y)
  }else if(isTRUE(subsample)){
    n_subsamples <- 1000
  }else{
    n_subsamples <- subsample
  }
  n_subsamples <- min(n_subsamples, length(y))

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
  if(max(c(y, -Inf)) < length(y) * 10){
     tab <- make_table(y)
  }else{
     tab <- list(numeric(0), numeric(0))
  }

  if(all(y == 0)){
    return_value$message <- "All counts y are 0."
    return_value$estimate <- 0
    return(return_value)
  }

  # Mu = 0 makes problems
  mean_vector[mean_vector == 0] <- 1e-6

  far_left_value <- conventional_score_function_fast(y, mu = mean_vector, log_theta = log(1e-8),
                                   model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment)
  if(far_left_value < 0){
    return_value$message <-  "Even for very small theta, no maximum identified"
    return_value$estimate <- 0
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
         })

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



#' Call gp_overdispersion_mle multiple times
#'
#' @return a vector with the overdispersion estimates per gene
#'
#' @keywords internal
estimate_overdispersions <- function(Y, mean_matrix, model_matrix, do_cox_reid_adjustment, subsample, verbose = FALSE){
  stopifnot(is.numeric(subsample))
  if(subsample < ncol(Y)){
    if(verbose){ message("Subsample data to ", subsample, " columns.") }
  }

  ## The estimate_overdispersions_fast() method is actually just doing the same
  ## as this vapply loop. However, the beachmat caching (?) is speeding up the procedure
  ## for HDF5 backed matrices by a factor of 50 and gives almost equivalent speed to in RAM
  ## methods.
  # vapply(seq_len(nrow(Y)), function(gene_idx){
  #   gampoi_overdispersion_mle(y = Y[gene_idx, ], mean_vector = mean_matrix[gene_idx, ],
  #                             model_matrix = model_matrix, do_cox_reid_adjustment = do_cox_reid_adjustment,
  #                             n_subsamples = n_subsamples)$estimate
  # }, FUN.VALUE = 0.0)
  estimate_overdispersions_fast(Y, mean_matrix, model_matrix, do_cox_reid_adjustment, subsample)

}




