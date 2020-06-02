
#' Make a quick first guess where reasonable beta would be
#'
#' @return a matrix with one column for each coefficient
#'
#' @keywords internal
estimate_betas_roughly <- function(Y, model_matrix, offset_matrix, pseudo_count = 1){
  if(nrow(Y) == 0) return(matrix(numeric(0), nrow = 0, ncol = ncol(model_matrix)))
  qrx <- qr(model_matrix)
  Q <- qr.Q(qrx)
  R <- qr.R(qrx)

  norm_log_count_mat <- t(log((Y / exp(offset_matrix) + pseudo_count)))
  t(solve(R, as.matrix(t(Q) %*% norm_log_count_mat)))
}


#' Estimate the Betas for Fixed Dispersions
#'
#' @return a list with two elements
#'   * `Beta` a matrix with one column for each coefficient
#'   * `iterations` the number of iterations
#'
#' @keywords internal
estimate_betas_fisher_scoring <- function(Y, model_matrix, offset_matrix,
                                          dispersions, beta_mat_init){
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))

  # The decision p > 30 is a rough heuristic:
  # For large p calculating the full info matrix (Xt W X) is not
  # worth the additional precision and it is faster to
  # approximate the info matrix with the diagonal element of Xt W X
  if(ncol(model_matrix) < 30){
    betaRes <- fitBeta_fisher_scoring(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init,
                                      ridge_penalty = 0, tolerance = 1e-8, max_iter =  1000)
  }else{
    # This one fails if there is considerable colinearity between columns
    betaRes <- fitBeta_diagonal_fisher_scoring(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init,
                                               tolerance = 1e-8, max_iter =  5000)
  }
  list(Beta = betaRes$beta_mat, iterations = betaRes$iter, deviances = betaRes$deviance)
}


#' Make a quick first guess where reasonable beta would be for an individual group
#'
#' @return a vector with the intercepts for each gene
#'
#' @keywords internal
estimate_betas_roughly_one_group <- function(Y, offset_matrix){
  log(DelayedMatrixStats::rowMeans2(Y / exp(offset_matrix)))
}


#' Estimate the Betas for Fixed Dispersions
#'
#' @return a list with two elements
#'   * `Beta` a matrix with one column and the intercept for each gene
#'   * `iterations` the number of iterations from the Newton-Raphson method
#'
#' @keywords internal
estimate_betas_one_group <- function(Y, offset_matrix,  dispersions, beta_vec_init){
  stopifnot(length(beta_vec_init) == nrow(Y))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))


  betaRes <- fitBeta_one_group(Y, offset_matrix, thetas = dispersions,
                               beta_start_values = beta_vec_init,
                               tolerance = 1e-8, maxIter = 100)

  list(Beta = matrix(betaRes$beta, nrow = nrow(Y), ncol = 1),
       iterations = betaRes$iter,
       deviances = betaRes$deviance)
}
