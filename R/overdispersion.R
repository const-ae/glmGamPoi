


#' Estimate a Single Overdispersion for a Vector of Counts
#'
#' @export
gp_overdispersion_mle <- function(y, mean_vector = mean(y),
                           model_matrix = matrix(1, nrow = length(y), ncol = 1),
                           do_cox_reid_adjustment = TRUE,
                           verbose = FALSE){

}




#' Call gp_overdispersion_mle multiple times
#'
#' Not exported
estimate_overdispersions <- function(Y, mean_matrix, model_matrix, do_cox_reid_adjustment){
  stopifnot(all(!missing(c(Y, mean_matrix, model_matrix, do_cox_reid_adjustment))))
}



