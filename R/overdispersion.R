


#' Estimate a Single Overdispersion for a Vector of Counts
#'
#' @export
gampoi_overdispersion_mle <- function(y, mean_vector = mean(y),
                           model_matrix = matrix(1, nrow = length(y), ncol = 1),
                           do_cox_reid_adjustment = TRUE,
                           verbose = FALSE){

  if(length(mean_vector) == 1){
    mean_vector <- rep(mean_vector, length(y))
  }
  stopifnot(length(y) == length(mean_vector))
  stopifnot(all(y >= 0))
  stopifnot(all(mean_vector >= 0))
  stopifnot(all(! is.na(y)))   # Cannot handle missing values
  stopifnot(all(! is.na(mean_vector)))
  stopifnot(all(is.finite(y)))
  stopifnot(all(is.finite(mean_vector)))

  # Decide if I use the Bandara approach or classical MLE
  # Rule is a rough heuristic
  if(max(y) < length(y)){
    # Do Bandara
    bandara_overdispersion_mle(y, mean_vector = mean_vector,
                               model_matrix = model_matrix,
                               do_cox_reid_adjustment = do_cox_reid_adjustment,
                               verbose = verbose)
  }else{
    # Do conventional optimization
    conventional_overdispersion_mle(y, mean_vector = mean_vector,
                                    model_matrix = model_matrix,
                                    do_cox_reid_adjustment = do_cox_reid_adjustment,
                                    verbose = verbose)
  }
}

bandara_overdispersion_mle <- function(y, mean_vector,
                           model_matrix = matrix(1, nrow = length(y), ncol = 1),
                           do_cox_reid_adjustment = TRUE,
                           verbose = FALSE){
  return_value = list(root = NA_real_, iterations = NA_real_, method = "bandara", message = "")


  if(all(y == 0)){
    return_value$message <- "All counts y are 0."
    return_value$root <- 0
    return(return_value)
  }

  # Mu = 0 makes problems
  mean_vector[mean_vector == 0] <- 1e-6

  # Common thing between all function calls
  # For explanation, see Bandara et al. (2019)
  cslv <- makeCumSumLookupVector(round(y))  # Cannot handle non integer values

  far_right_value <- score_function_bandara_fast(y, cumsumLookupTable = cslv,
                                                 mean_vector, r = 1000,
                                                 model_matrix = model_matrix,
                                                 do_cr_adj = do_cox_reid_adjustment)
  if(far_right_value > 0){
    return_value$message <- "Even for very small theta, no maximum identified"
    return_value$root <- 0
    return(return_value)
  }

  my <- mean(y)
  bvar <- sum((y-my)^2) / length(y)
  rmme <- my / (bvar - my)
  start_pos <- NA
  gamma_factor <- 0.9
  if(is.na(rmme) || (rmme < 0 && max(y) > 1)){
    # Exceptional case
    N1 <- sum(y >= 1)
    start_13 <- (-3 + sqrt(24 * length(y) * my / N1 - 15)) / (2 * N1 / (length(y) * my))
    start_pos <- start_13 / gamma_factor
  }

  if(rmme < 0){
    rmme <- 1e-3
  }
  iter <- 1
  repeat{
    start_pos <- rmme * gamma_factor
    if(score_function_bandara_fast(y, cumsumLookupTable = cslv,
                                   mean_vector, r= start_pos,
                                   model_matrix = model_matrix,
                                   do_cr_adj = do_cox_reid_adjustment) > 0){
      break
    }else if(iter >= 100){
      return_value$message <- paste0("Couldn't find starting position. Stopped after ", iter, " iterations.")
      return(return_value)
    }else{
      gamma_factor <- gamma_factor / 2
      iter <- iter + 1
    }
  }

  root_info <- pracma::newtonRaphson(function(r){
    score_function_bandara_fast(y, cumsumLookupTable = cslv,
                                mu = mean_vector, r = r,
                                model_matrix = model_matrix,
                                do_cr_adj = do_cox_reid_adjustment)
  }, x0 = start_pos,
  dfun = function(r){
    score_deriv_function_bandara_fast(y, cumsumLookupTable = cslv,
                                      mu = mean_vector, r= r,
                                      model_matrix = model_matrix,
                                      do_cr_adj = do_cox_reid_adjustment)
  }, tol = .Machine$double.eps^0.25)

  return_value$root <- 1/root_info$root
  return_value$iterations <- root_info$niter
  return_value$message <- "success"
  return_value
}


conventional_overdispersion_mle <- function(y, mean_vector,
                                       model_matrix = matrix(1, nrow = length(y), ncol = 1),
                                       do_cox_reid_adjustment = TRUE,
                                       verbose = FALSE){
  return_value = list(root = NA_real_, iterations = NA_real_, method = "conventional", message = "")


  if(all(y == 0)){
    return_value$message <- "All counts y are 0."
    return_value$root <- 0
    return(return_value)
  }

  # Mu = 0 makes problems
  mean_vector[mean_vector == 0] <- 1e-6

  far_left_value <- conventional_score_function_fast(y, mu = mean_vector, log_theta = log(1/1000),
                                   model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment)
  if(far_left_value < 0){
    return_value$message <-  "Even for very small theta, no maximum identified"
    return_value$root <- 0
    return(return_value)
  }

  mu <- mean(y)
  start_value <- (var(y) - mu) / mu^2
  if(is.na(start_value) || start_value < 0){
    start_value <- 0.5
  }

  res <- nlminb(start = log(start_value),
         objective = function(log_theta){
           - conventional_loglikelihood_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment)
         }, gradient = function(log_theta){
           - conventional_score_function_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment)
         }, hessian = function(log_theta){
           res <- conventional_deriv_score_function_fast(y, mu = mean_vector, log_theta = log_theta,
                             model_matrix = model_matrix, do_cr_adj = do_cox_reid_adjustment)
           matrix(- res, nrow = 1, ncol = 1)
         })

  return_value$root <- exp(res$par)
  return_value$iterations <- res$iterations
  return_value$message <- res$message
  return_value
}



#' Call gp_overdispersion_mle multiple times
#'
#' Not exported
#' @keywords internal
estimate_overdispersions <- function(Y, mean_matrix, model_matrix, do_cox_reid_adjustment){

  vapply(seq_len(nrow(Y)), function(gene_idx){
    gampoi_overdispersion_mle(y = Y[gene_idx, ], mean_vector = mean_matrix[gene_idx, ],
                              model_matrix = model_matrix, do_cox_reid_adjustment = do_cox_reid_adjustment)$root
  }, FUN.VALUE = 0.0)

}




estimate_dispersion_by_moment <- function(Y, model_matrix, offset_matrix){
  xim <- 1/mean(matrixStats::colMeans2(exp(offset_matrix)))
  bv <- matrixStats::rowVars(Y)
  bm <- matrixStats::rowMeans2(Y)
  (bv - xim * bm) / bm^2
}

estimate_dispersions_roughly <- function(Y, model_matrix, offset_matrix){
  roughDisp <- DESeq2:::roughDispEstimate(y = Y / exp(offset_matrix),
                                          x = model_matrix)
  momentsDisp <- estimate_dispersion_by_moment(Y, model_matrix, offset_matrix)
  disp_rough <- pmin(roughDisp, momentsDisp)
  ifelse(is.na(disp_rough) | disp_rough < 0, 0, disp_rough)
}

