
#' Make a quick first guess where reasonable beta would be
#'
#' @return a matrix with one column for each coefficient
#'
#' @keywords internal
estimate_betas_roughly <- function(Y, model_matrix, offset_matrix, pseudo_count = 1, ridge_penalty = NULL){
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))
  if(is.vector(offset_matrix, mode = "numeric")){
    stopifnot(length(offset_matrix) == ncol(Y))
  }else{
    stopifnot(dim(offset_matrix) == dim(Y))
  }

  if(nrow(Y) == 0){
    return(matrix(numeric(0), nrow = 0, ncol = ncol(model_matrix)))
  }

  if(is.null(ridge_penalty)){
    qrx <- qr(model_matrix)
  }else if(is.matrix(ridge_penalty)){
    qrx <- qr(rbind(model_matrix, ridge_penalty))
  }else if(is.numeric(ridge_penalty)){
    qrx <- qr(rbind(model_matrix, diag(ridge_penalty, nrow = length(ridge_penalty))))
  }else{
    stop("Illegal ridge penalty definition")
  }


  Q <- qr.Q(qrx)[seq_len(nrow(model_matrix)),,drop=FALSE]
  R <- qr.R(qrx)

  if(is.vector(offset_matrix, mode = "numeric")){
    norm_Y <- div_mtx_colwise(Y, exp(offset_matrix))
  }else{
    norm_Y <- div_mtx_elemwise(Y, exp(offset_matrix))
  }

  if (pseudo_count == 1) {
    norm_log_count_mat <- log1p(norm_Y)
  } else {
    norm_log_count_mat <- log(norm_Y + pseudo_count)
  }

  t(solve(R, as.matrix(Matrix::tcrossprod(t(Q), norm_log_count_mat))))
}


#' Estimate the Betas for Fixed Dispersions
#'
#' @return a list with two elements
#'   * `Beta` a matrix with one column for each coefficient
#'   * `iterations` the number of iterations
#'
#' @keywords internal
#' @importFrom beachmat initializeCpp
estimate_betas_fisher_scoring <- function(Y, model_matrix, offset_matrix,
                                          dispersions, beta_mat_init, ridge_penalty,
                                          try_recovering_convergence_problems = TRUE,
                                          max_iter = 1000,
                                          do_parallel = 0){
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  if(is.vector(offset_matrix, mode = "numeric")){
    stopifnot(length(offset_matrix) == ncol(Y))
  }else{
    stopifnot(dim(offset_matrix) == dim(Y))
  }

  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))

  if(! is.null(ridge_penalty) && ! is.matrix(ridge_penalty)){
    ridge_target <- attr(ridge_penalty, "target")
    ridge_penalty <- diag(ridge_penalty, nrow = length(ridge_penalty))
    attr(ridge_penalty, "target") <- ridge_target
  }

  exp_offset_matrix <- exp(offset_matrix)
  if(is.vector(exp_offset_matrix, mode = "numeric")){
    exp_offset_matrix <- matrix(exp_offset_matrix, nrow = 1)
  }
  betaRes <- fitBeta_fisher_scoring(initializeCpp(Y), model_matrix, initializeCpp(exp_offset_matrix), dispersions, beta_mat_init,
                                    ridge_penalty_nl = ridge_penalty, tolerance = 1e-8,
                                    max_rel_mu_change = 1e5, max_iter = max_iter, try_recov_w_optim = try_recovering_convergence_problems, do_parallel = do_parallel)
  warn_non_convergence(betaRes$iter == max_iter, rownames(Y))

  list(Beta = betaRes$beta_mat, iterations = betaRes$iter, deviances = betaRes$deviance)
}

warn_non_convergence <- function(not_converged, rownames){
  if(any(not_converged)){
    # Estimate didn't converge for some gene :(
    labels <- if(! is.null(rownames)){
      rownames[not_converged]
    }else{
      which(not_converged)
    }
    warning("Beta estimation did not converge for ", paste0(head(labels), collapse = ", "),
            if(length(labels) > 6){", ..."}, ".\n",
            "Will continue anyways and ignore those rows in subsequent calls.")
  }
}

estimate_betas_optim <- function(Y, model_matrix, offset_matrix, dispersions, beta_mat_init, ridge_penalty, max_iter = 1000, do_parallel = 0){
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  if(is.vector(offset_matrix, mode = "numeric")){
    stopifnot(length(offset_matrix) == ncol(Y))
  }else{
    stopifnot(dim(offset_matrix) == dim(Y))
  }
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))


  if(! is.null(ridge_penalty) && ! is.matrix(ridge_penalty)){
    ridge_target <- attr(ridge_penalty, "target")
    ridge_penalty <- diag(ridge_penalty, nrow = length(ridge_penalty))
    attr(ridge_penalty, "target") <- ridge_target
  }

  fitBeta_optim(initializeCpp(Y), model_matrix, initializeCpp(exp(offset_matrix)), dispersions, beta_mat_init, ridge_penalty, max_iter, do_parallel = do_parallel)
}


#' Make a quick first guess where reasonable beta would be for a set of groups
#'
#' @return a matrix with the mean per group for each gene
#'
#' @keywords internal
estimate_betas_roughly_group_wise <- function(Y, offset_matrix, groups){
  if(is.vector(offset_matrix, mode = "numeric")){
    norm_Y <- div_mtx_colwise(Y, exp(offset_matrix))
  } else {
    norm_Y <- div_mtx_elemwise(Y, exp(offset_matrix))
  }
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
#' @importFrom beachmat initializeCpp
estimate_betas_group_wise <- function(Y, offset_matrix,  dispersions, beta_group_init = NULL, beta_mat_init = NULL, groups, model_matrix, max_iter = 100, do_parallel = 0){
  stopifnot(nrow(beta_group_init) == nrow(Y))
  stopifnot(ncol(beta_group_init) == length(unique(groups)))
  stopifnot(length(dispersions) == nrow(Y))
  if(is.vector(offset_matrix, mode = "numeric")){
    stopifnot(length(offset_matrix) == ncol(Y))
  }else{
    stopifnot(dim(offset_matrix) == dim(Y))
  }
  stopifnot(is.null(beta_mat_init) != is.null(beta_group_init))
  if(is.null(beta_group_init)){
    # Calculate group_init based on Beta
    first_occurence_in_groups <- match(unique(groups), groups)
    beta_group_init <- Matrix::tcrossprod(beta_mat_init, model_matrix[first_occurence_in_groups, ,drop=FALSE])
  }

  Beta_res_list <- lapply(unique(groups), function(gr){
    chosen <- gr == groups
    Y_gr <- Y[, chosen, drop = FALSE]
    if(is.vector(offset_matrix, mode = "numeric")){
      offset_gr <- matrix(offset_matrix[chosen, drop = FALSE], nrow = 1)
    }else{
      offset_gr <- offset_matrix[, chosen, drop = FALSE]
    }
    fitBeta_one_group(initializeCpp(Y_gr),
                                 initializeCpp(offset_gr), thetas = dispersions,
                                 beta_start_values = beta_group_init[, gr == unique(groups),drop=TRUE],
                                 tolerance = 1e-8, max_iter = max_iter, do_parallel = do_parallel)
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
       deviances = matrixStats::rowSums2(Deviance_mat))
}

estimate_betas_group_wise_optimize_helper <- function(y, offset, theta, lower_bound = -30, upper_bound = 30){
  optimize(function(beta){
    sum(dnbinom(y, mu = exp(beta + offset), size = 1/theta, log = TRUE))
  }, lower = lower_bound, upper = upper_bound, maximum = TRUE)$maximum
}

