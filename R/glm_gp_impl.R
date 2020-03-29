


#' Internal Function to Fit a Gamma-Poisson GLM
#'
#' @inheritParams glm_gp
#' @inheritParams gampoi_overdispersion_mle
#' @param Y any matrix-like object (e.g. `matrix()`, `DelayedArray()`, `HDF5Matrix()`) with
#'   one column per sample and row per gene.
#'
#' @return a list with four elements
#'  * `Beta` the coefficient matrix
#'  * `overdispersion` the vector with the estimated overdispersions
#'  * `Mu` a matrix with the corresponding means for each gene
#'     and sample
#'  * `size_factors` a vector with the size factor for each
#'    sample
#'
#' @seealso [glm_gp()] and [gampoi_overdispersion_mle()]
#' @keywords internal
glm_gp_impl <- function(Y, model_matrix,
                        offset = 0,
                        size_factors = TRUE,
                        overdispersion = TRUE,
                        do_cox_reid_adjustment = TRUE,
                        subsample = FALSE,
                        verbose = FALSE){
  if(is.vector(Y)){
    Y <- matrix(Y, nrow = 1)
  }
  # Error conditions
  stopifnot(is.matrix(Y) || is(Y, "DelayedArray"))
  stopifnot(is.matrix(model_matrix) && nrow(model_matrix) == ncol(Y))
  validate_Y_matrix(Y)
  subsample <- handle_subsample_parameter(Y, subsample)

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
    Beta <- estimate_betas_one_group(Y, offset_matrix = offset_matrix,
                                         dispersions = disp_init, beta_vec_init = beta_vec_init)$Beta
  }else{
    # Init beta with reasonable values
    if(verbose){ message("Make initial beta estimate") }
    beta_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix)
    if(verbose){ message("Estimate beta") }
    Beta <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                              dispersions = disp_init, beta_mat_init = beta_init)$Beta
  }

  # Calculate corresponding predictions
  # Mu <- exp(Beta %*% t(model_matrix) + offset_matrix)
  Mu <- calculate_mu(Beta, model_matrix, offset_matrix)

  # Make estimate of over-disperion
  if(isTRUE(overdispersion)){
    if(verbose){ message("Estimate dispersion") }
    disp_est <- estimate_overdispersions(Y, Mu, model_matrix = model_matrix,
                                         do_cox_reid_adjustment = do_cox_reid_adjustment,
                                         subsample = subsample, verbose = verbose)

    # Estimate the betas again (only necessary if disp_est has changed)
    if(verbose){ message("Estimate beta again") }
    if(only_intercept_model){
      Beta <- estimate_betas_one_group(Y, offset_matrix = offset_matrix,
                                       dispersions = disp_est, beta_vec_init = Beta[,1])$Beta
    }else{
      Beta <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                            dispersions = disp_est, beta_mat_init = Beta)$Beta
    }

    # Calculate corresponding predictions
    Mu <- calculate_mu(Beta, model_matrix, offset_matrix)

  }else{
    # Use disp_init, because it is already in vector shape
    disp_est <- disp_init
  }


  # Return everything
  list(Beta = Beta, overdispersions = disp_est,
       Mu = Mu, size_factors = size_factors)
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


