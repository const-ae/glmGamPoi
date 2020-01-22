


#' Low-level Function to Fit a Gamma-Poisson GLM
#'
#' @export
glm_gp_impl <- function(Y, design_matrix,
                        offset = 0,
                        size_factors = TRUE,
                        overdispersion = NULL,
                        do_cox_reid_adjustment = TRUE,
                        verbose = FALSE){
  # Error conditions
  stopifnot(is.matrix(Y))
  stopifnot(is.matrix(design_matrix) && nrow(design_matrix) == ncol(Y))

  # Combine offset and size factor
  off_and_sf <- combine_size_factors_and_offset(offset, size_factors, Y, verbose = verbose)
  offset_matrix <- off_and_sf$offset_matrix
  size_factors <- off_and_sf$size_factors

  # Decide if there is only the intercept
  only_intercept_model <- ncol(design_matrix) == 1 && all(design_matrix == 1)

  # If no overdispersion, make rough first estimate
  if(is.null(overdispersion)){
    if(verbose){ message("Make initial dispersion estimate") }
    disp_init <- estimate_dispersions_roughly(Y, design_matrix, offset_matrix = offset_matrix)
  }else{
    stopifnot(is.numeric(overdispersion) && (length(overdispersion) == 1 || length(overdispersion) == nrow(Y)))
    if(length(overdispersion) == 1){
      disp_init <- rep(overdispersion, times = nrow(Y))
    }else{
      disp_init <- overdispersion
    }
  }


  # Estimate the betas
  if(verbose){ message("Estimate beta") }
  if(only_intercept_model){
    Beta_est <- estimate_betas_one_group(Y, offset_matrix = offset_matrix, dispersion = disp_init)$Beta
  }else{
    # Init beta with reasonable values
    if(verbose){ message("Make rough initial beta estimate") }
    beta_init <- estimate_betas_roughly(Y, design_matrix, offset_matrix = offset_matrix)
    Beta_est <- estimate_betas(Y, model_matrix = design_matrix, offset_matrix = offset_matrix,
                               dispersions = disp_init, beta_mat_init = beta_init)$Beta
  }

  # Calculate corresponding predictions
  Mu_est <- exp(Beta_est %*% t(design_matrix) + offset_matrix)

  # Make estimate of over-disperion
  if(is.null(overdispersion)){
    if(verbose){ message("Estimate dispersion") }
    disp_est <- estimate_overdispersions(Y, Mu_est, model_matrix = design_matrix,
                                         do_cox_reid_adjustment = TRUE)
  }else{
    disp_est <- disp_init
  }

  # Estimate the betas again
  if(verbose){ message("Estimate beta again") }
  if(only_intercept_model){
    Beta_est <- estimate_betas_one_group(Y, offset_matrix = offset_matrix,
                                         dispersion = disp_est, beta_vec_init = Beta_est[,1])$Beta
  }else{
    Beta_est <- estimate_betas(Y, model_matrix = design_matrix, offset_matrix = offset_matrix,
                               dispersion = disp_est, beta_mat_init = Beta_est)$Beta
  }

  # Calculate corresponding predictions
  mu_mat <- exp(Beta_est %*% t(design_matrix) + offset_matrix)

  # Return everything
  list(Beta_est = Beta_est, overdispersions = disp_est,
       Mu_est = Mu_est, size_factors = size_factors)
}


combine_size_factors_and_offset <- function(offset, size_factors, Y, verbose = FALSE){
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)
  if(is.matrix(offset)){
    stopifnot(dim(offset) == c(n_genes, n_samples))
    offset_matrix <- offset
  }else{
    stopifnot(length(offset) == 1 || length(offset) == n_samples)
    offset_matrix <- matrix(offset, nrow=n_genes, ncol = n_samples, byrow = TRUE)
  }
  if(isTRUE(size_factors)){
    if(verbose){ message("Calculate Size Factors") }
    lsf <- log(estimate_size_factors(Y))
  }else if(isFALSE(size_factors)){
    lsf <- 0
  }else{
    stopifnot(is.numeric(size_factors) && (length(size_factors) == 1 || length(size_factors) == n_samples))
    if(any(size_factors < 0)){
      stop("size factor 'size_factors' must be larger than 0")
    }
    if(length(size_factors) == 1){
      lsf <- rep(log(size_factors), n_samples)
    }else{
      lsf <- log(size_factors)
    }
  }
  lsf_mat <- matrix(lsf, nrow = n_genes, ncol = n_samples, byrow = TRUE)
  offset_matrix <- offset_matrix + lsf_mat
  list(offset_matrix = offset_matrix, size_factors = exp(lsf))
}

