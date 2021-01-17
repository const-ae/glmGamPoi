
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
                                          dispersions, beta_mat_init, ridge_penalty){
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && nrow(ridge_penalty) == ncol(model_matrix) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))

  if(! is.null(ridge_penalty) && ! is.matrix(ridge_penalty)){
    ridge_penalty <- diag(ridge_penalty, nrow = length(ridge_penalty))
  }

  betaRes <- fitBeta_fisher_scoring(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init,
                                    ridge_penalty = ridge_penalty, tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  1000)

  list(Beta = betaRes$beta_mat, iterations = betaRes$iter, deviances = betaRes$deviance)
}


#' Make a quick first guess where reasonable beta would be for a set of groups
#'
#' @return a matrix with the mean per group for each gene
#'
#' @keywords internal
estimate_betas_roughly_group_wise <- function(Y, offset_matrix, groups){
  norm_Y <- Y / exp(offset_matrix)
  do.call(cbind, lapply(unique(groups), function(gr){
    log(DelayedMatrixStats::rowMeans2(norm_Y, cols = groups == gr))
  }))
}


#' Estimate the Betas for Fixed Dispersions
#'
#' @return a list with three elements
#'   * `Beta` a matrix with one column per group and a row for each gene
#'   * `iterations` the number of iterations from the Newton-Raphson method
#'   * `deviances` the deviance for each gene (sum of the deviance per group)
#'
#' @keywords internal
estimate_betas_group_wise <- function(Y, offset_matrix,  dispersions, beta_group_init = NULL, beta_mat_init = NULL, groups, model_matrix){
  stopifnot(nrow(beta_group_init) == nrow(Y))
  stopifnot(ncol(beta_group_init) == length(unique(groups)))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))
  stopifnot(is.null(beta_mat_init) != is.null(beta_group_init))
  if(is.null(beta_group_init)){
    # Calculate group_init based on Beta
    first_occurence_in_groups <- match(unique(groups), groups)
    beta_group_init <- beta_mat_init %*% t(model_matrix[first_occurence_in_groups, ,drop=FALSE])
  }

  Beta_res_list <- lapply(unique(groups), function(gr){
    betaRes <- fitBeta_one_group(Y[, gr == groups, drop = FALSE],
                                 offset_matrix[, gr == groups, drop = FALSE], thetas = dispersions,
                                 beta_start_values = beta_group_init[, gr == unique(groups),drop=TRUE],
                                 tolerance = 1e-8, maxIter = 100)
  })
  Beta <- do.call(cbind, lapply(Beta_res_list, function(x) x$beta))
  Iteration_mat <- do.call(cbind, lapply(Beta_res_list, function(x) x$iter))
  Deviance_mat <- do.call(cbind, lapply(Beta_res_list, function(x) x$deviance))

  # How about rotating the Beta into the right place?!
  Beta <- pmax(Beta, -1e8)
  first_occurence_in_groups <- match(unique(groups), groups)
  if(nrow(Beta) > 0){
    Beta <- t(solve(model_matrix[first_occurence_in_groups, ,drop=FALSE], t(Beta)))
  }

  list(Beta = Beta,
       iterations = matrixStats::rowSums2(Iteration_mat),
       deviances = matrixStats::rowWeightedMeans(Deviance_mat, w = c(table(groups))))
}

estimate_betas_group_wise_optimize_helper <- function(y, offset, theta, lower_bound = -30, upper_bound = 30){
  optimize(function(beta){
    sum(dnbinom(y, mu = exp(beta + offset), size = 1/theta, log = TRUE))
  }, lower = lower_bound, upper = upper_bound, maximum = TRUE)$maximum
}

