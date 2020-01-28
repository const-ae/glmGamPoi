


#' Low-level Function to Fit a Gamma-Poisson GLM
#'
#' @export
glm_gp_impl <- function(Y, model_matrix,
                        offset = 0,
                        size_factors = TRUE,
                        overdispersion = TRUE,
                        do_cox_reid_adjustment = TRUE,
                        n_subsamples = min(1000, ncol(Y)),
                        verbose = FALSE){
  if(is.vector(Y)){
    Y <- matrix(Y, nrow = 1)
  }
  # Error conditions
  stopifnot(is.matrix(Y) || is(Y, "DelayedArray"))
  stopifnot(is.matrix(model_matrix) && nrow(model_matrix) == ncol(Y))
  validate_Y_matrix(Y)


  # Combine offset and size factor
  off_and_sf <- combine_size_factors_and_offset(offset, size_factors, Y, verbose = verbose)
  offset_matrix <- off_and_sf$offset_matrix
  size_factors <- off_and_sf$size_factors

  # Decide if there is only the intercept
  only_intercept_model <- ncol(model_matrix) == 1 && all(model_matrix == 1)

  # If no overdispersion, make rough first estimate
  if(isTRUE(overdispersion)){
    if(verbose){ message("Make initial dispersion estimate") }
    disp_init <- estimate_dispersions_roughly(Y, model_matrix, offset_matrix = offset_matrix)
  }else if(isFALSE(overdispersion)){
    disp_init <- rep(0, times = nrow(Y))
  }else{
    stopifnot(is.numeric(overdispersion) && (length(overdispersion) == 1 || length(overdispersion) == nrow(Y)))
    if(length(overdispersion) == 1){
      disp_init <- rep(overdispersion, times = nrow(Y))
    }else{
      disp_init <- overdispersion
    }
  }


  # Estimate the betas
  if(only_intercept_model){
    if(verbose){ message("Make initial beta estimate") }
    beta_vec_init <- estimate_betas_roughly_one_group(Y, offset_matrix)
    if(verbose){ message("Estimate beta") }
    Beta_est <- estimate_betas_one_group(Y, offset_matrix = offset_matrix,
                                         dispersions = disp_init, beta_vec_init = beta_vec_init)$Beta
  }else{
    # Init beta with reasonable values
    if(verbose){ message("Make initial beta estimate") }
    beta_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix)
    if(verbose){ message("Estimate beta") }
    Beta_est <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                              dispersions = disp_init, beta_mat_init = beta_init)$Beta
  }

  # Calculate corresponding predictions
  # Mu_est <- exp(Beta_est %*% t(model_matrix) + offset_matrix)
  Mu_est <- calculate_mu(Beta_est, model_matrix, offset_matrix)

  # Make estimate of over-disperion
  if(isTRUE(overdispersion)){
    if(verbose){ message("Estimate dispersion") }
    disp_est <- estimate_overdispersions(Y, Mu_est, model_matrix = model_matrix,
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
                                         dispersions = disp_est, beta_vec_init = Beta_est[,1])$Beta
  }else{
    Beta_est <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                              dispersions = disp_est, beta_mat_init = Beta_est)$Beta
  }

  # Calculate corresponding predictions
  Mu_est <- calculate_mu(Beta_est, model_matrix, offset_matrix)

  # Return everything
  list(Beta_est = Beta_est, overdispersions = disp_est,
       Mu_est = Mu_est, size_factors = size_factors)
}



validate_Y_matrix <- function(Y){
  which_neg_values <- which(DelayedMatrixStats::colAnys(Y < 0))
  if(length(which_neg_values)){
    stop("The data contains negative values (Y < 0) in columns: ",
         if(length(which_neg_values) > 5) paste0(paste(head(which_neg_values), collapse = ","), ", ...")
         else paste(head(which_neg_values), collapse = ","), "\n",
         "This is does not make sense as input data which has a support from (0 to 2^31-1).\n", call. = FALSE)
  }
  cs <- DelayedMatrixStats::colSums2(Y)
  which_is_na <- which(is.na(cs))
  if(length(which_is_na)){
    stop("The data contains missing values ('NA') in columns: ",
           if(length(which_is_na) > 5) paste0(paste(head(which_is_na), collapse = ","), ", ...")
           else paste(head(which_is_na), collapse = ","), "\n",
         "This is currently not supported by glmGamPoi.\n", call. = FALSE)
  }
  which_is_inf <- which(is.infinite(cs))
  if(length(which_is_inf)){
    stop("The data contains infinite values ('Inf' or '-Inf') in columns: ",
         if(length(which_is_inf) > 5) paste0(paste(head(which_is_inf), collapse = ","), ", ...")
         else paste(head(which_is_inf), collapse = ","), "\n",
         "This is currently not supported by glmGamPoi.\n", call. = FALSE)
  }
}


