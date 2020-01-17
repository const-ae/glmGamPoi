


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

  # Decide if I use the Bandara approach or classical MLE
  if(max(y) < length(y)){
    # Do Bandara
    bandara_overdispersion_mle(y, mean_vector = mean_vector,
                               model_matrix = model_matrix,
                               do_cox_reid_adjustment = do_cox_reid_adjustment,
                               verbose = verbose)
  }else{
    # Do conventional optimization
    conventional_overdispersion_mle(y, mean_vector = mean_vecot,
                                    model_matrix = model_matrix,
                                    do_cox_reid_adjustment = do_cox_reid_adjustment,
                                    verbose = verbose)
  }
}

bandara_overdispersion_mle <- function(y, mean_vector,
                           model_matrix = NULL,
                           do_cox_reid_adjustment = TRUE,
                           verbose = FALSE){
  return_value = list(root = NA_real_, iterations = NA_real_, method = "bandara", message = "")

  # Common thing between all function calls
  # For explanation, see Bandara et al. (2019)
  cslv <- makeCumSumLookupVector(y)

  far_right_value <- score_function_bandara_fast(y, cumsumLookupTable = cslv,
                                                 mean_vector, r = 1000,
                                                 model_matrix = model_matrix,
                                                 do_cr_adj = do_cox_reid_adjustment)
  if(far_right_value > 0){
    return_value$message <- "Even for very large theta, no maximum identified"
    return(return_value)
  }

  my <- mean(y)
  bvar <- sum((y-my)^2) / length(y)
  rmme <- my / (bvar - my)
  start_pos <- NA
  gamma_factor <- 0.9
  if(rmme < 0 && max(y) > 1){
    # Exceptional case
    N1 <- sum(y >= 1)
    start_13 <- (-3 + sqrt(24 * length(y) * mean(y) / N1 - 15)) / (2 * N1 / (length(y) * mean(y)))
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
                                      mu = mean_vector, theta=1/r,
                                      model_matrix = model_matrix,
                                      do_cr_adj = do_cox_reid_adjustment)
  }, tol = .Machine$double.eps^0.25)
  return_value$root <- 1/root_info$root
  return_value$iterations <- root_info$iter
  return_value$message <- "success"
  return_value
}


conventional_overdispersion_mle <- function(y, mean_vector,
                                       model_matrix = NULL,
                                       do_cox_reid_adjustment = TRUE,
                                       verbose = FALSE){
  return_value = list(root = NA_real_, iterations = NA_real_, method = "conventional", message = "")

  mu <- mean(y)
  start_value <- (var(y) - mu) / mu^2
  if(start_value < 0){
    return_value$message <- "Failed to find start value, mean too large"
    return(return_value)
  }

  nb_ll <- function(y, mu, theta, model_matrix, do_cox_reid_adjustment){
    if(do_cox_reid_adjustment){
      b <- t(model_matrix) %*% diag(1/(1/mu + theta), nrow = length(y)) %*% model_matrix
      cr_term <- 0.5 * log(det(b))
    }else{
      cr_term <- 0
    }
    ll_term <- sum(lgamma(y + 1/theta) - lgamma(1/theta) - y * log(mu + 1/theta) - 1/theta * log(1 + mu * theta))
    cr_term + ll_term
  }

  res <- nlminb(start = start_value,
         objective = function(theta){
           - nb_ll(y, mean_vector, exp(theta), model_matrix, do_cox_reid_adjustment)
         })
  return_value$root <- exp(res$par)
  return_value$iterations <- res$iterations
  return_value$message <- res$message
  return_value
}



#' Call gp_overdispersion_mle multiple times
#'
#' Not exported
estimate_overdispersions <- function(Y, mean_matrix, model_matrix, do_cox_reid_adjustment){
  stopifnot(all(!missing(c(Y, mean_matrix, model_matrix, do_cox_reid_adjustment))))
}



