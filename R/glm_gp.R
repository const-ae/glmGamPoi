

#' Fit a Gamma-Poisson Generalized Linear Model
#'
#' @export
glm_gp <- function(data, design = ~ 1,
                   col_data = NULL,
                   reference_level = NULL,
                   size_factors = TRUE,
                   overdispersion = TRUE,
                   do_cox_reid_adjustment = TRUE,
                   n_subsamples = min(1000, ncol(Y)),
                   verbose = FALSE){
  stop("Not yet implemented")
  # Convert the formula to a design_matrix
  # (pay attention to the names of design)

  # Check if the design_matrix is valid

  # Check if data is valid

  # Call glm_gp_impl()

  # Make sure that the output is nice and
  # beautiful

}




