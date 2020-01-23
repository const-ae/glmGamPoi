
#' Make a quick first guess where reasonable beta would be
#'
#' @importFrom Matrix t
#'
#' @keywords internal
estimate_betas_roughly <- function(Y, model_matrix, offset_matrix, pseudo_count = 1){
  qrx <- qr(model_matrix)
  Q <- qr.Q(qrx)
  R <- qr.R(qrx)

  norm_log_count_mat <- t(log((Y / exp(offset_matrix) + pseudo_count)))
  t(solve(R, as.matrix(t(Q) %*% norm_log_count_mat)))
}


#' Estimate the Betas for Fixed Dispersions
#'
#'
#' @keywords internal
estimate_betas <- function(Y, model_matrix, offset_matrix,
                           dispersions = NULL, beta_mat_init = NULL){

  if(is.null(beta_mat_init)){
    beta_mat_init<- estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix)
  }
  if(is.null(dispersions)){
    mu_mat <- exp(beta_mat_init %*% t(model_matrix) + offset_matrix)
    disps <- estimate_overdispersions(Y, mu_mat, model_matrix, do_cox_reid_adjustment = TRUE)
  }


  betaRes <- fitBeta(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init,
          tolSEXP = 1e-8, maxitSEXP =  100, minmuSEXP = 1e-6)

  list(Beta = betaRes$beta_mat, iterations = betaRes$iter)
}



#' Estimate the Betas for Fixed Dispersions
#'
#'
#' @keywords internal
estimate_betas_one_group <- function(Y, offset_matrix,  dispersion, beta_vec_init = NULL){
  if(is.null(beta_vec_init)){
    beta_vec_init <- log(rowMeans(Y / exp(offset_matrix)))
  }

  betaRes <- fitBeta_one_group(Y, offset_matrix, thetas = dispersion,
                               beta_start_values = beta_vec_init,
                               tolerance = 1e-8, maxIter = 100)


  list(Beta = matrix(betaRes$beta, nrow = nrow(Y), ncol = 1),
       iterations = betaRes$iter)
}
