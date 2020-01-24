


#' Low-level Function to Fit a Gamma-Poisson GLM
#'
#' @export
glm_gp_impl <- function(Y, design_matrix,
                        offset = 0,
                        size_factors = TRUE,
                        overdispersion = NULL,
                        do_cox_reid_adjustment = TRUE,
                        n_subsamples = min(1000, ncol(Y)),
                        verbose = FALSE){
  # Error conditions
  # stopifnot(is.matrix(Y))
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
  # Mu_est <- exp(Beta_est %*% t(design_matrix) + offset_matrix)
  Mu_est <- calculate_mu(Beta_est, design_matrix, offset_matrix)

  # Make estimate of over-disperion
  if(is.null(overdispersion)){
    if(verbose){ message("Estimate dispersion") }
    disp_est <- estimate_overdispersions(Y, Mu_est, model_matrix = design_matrix,
                                         do_cox_reid_adjustment = TRUE,
                                         n_subsamples = n_subsamples, verbose = verbose)
  }else{
    # Use disp_init, because it is already in vector shape
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
  Mu_est <- calculate_mu(Beta_est, design_matrix, offset_matrix)

  # Return everything
  list(Beta_est = Beta_est, overdispersions = disp_est,
       Mu_est = Mu_est, size_factors = size_factors)
}



