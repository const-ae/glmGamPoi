
#' Make a quick first guess where reasonable beta would be
#'
#' @return a matrix with one column for each coefficient
#'
#' @keywords internal
estimate_betas_roughly <- function(Y, model_matrix, offset_matrix, pseudo_count = 1, ridge_penalty = NULL){
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))

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
                                          dispersions, beta_mat_init, ridge_penalty,
                                          try_recovering_convergence_problems = TRUE){
  max_iter <- 1000
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))

  if(! is.null(ridge_penalty) && ! is.matrix(ridge_penalty)){
    ridge_target <- attr(ridge_penalty, "target")
    ridge_penalty <- diag(ridge_penalty, nrow = length(ridge_penalty))
    attr(ridge_penalty, "target") <- ridge_target
  }

  betaRes <- fitBeta_fisher_scoring(Y, model_matrix, exp(offset_matrix), dispersions, beta_mat_init,
                                    ridge_penalty_nl = ridge_penalty, tolerance = 1e-8,
                                    max_rel_mu_change = 1e5, max_iter =  max_iter)
  not_converged <- betaRes$iter == max_iter
  if(try_recovering_convergence_problems & any(not_converged)){
    # Try again with optim
    betaRes2 <- estimate_betas_optim(Y[not_converged,,drop=FALSE], model_matrix,
                                     offset_matrix[not_converged,,drop=FALSE],
                                     dispersions = dispersions[not_converged],
                                     beta_mat_init = beta_mat_init[not_converged,,drop=FALSE],
                                     ridge_penalty = ridge_penalty, max_iter = max_iter)
    betaRes$beta_mat[not_converged, ] <- betaRes2$Beta
    betaRes$deviance[not_converged] <- betaRes2$deviances
    betaRes$iter[not_converged] <- betaRes2$iterations
  }
  # Don't use 'not_converged' because optim might recover some cases
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

estimate_betas_optim <- function(Y, model_matrix, offset_matrix, dispersions, beta_mat_init, ridge_penalty, max_iter = 1000){
  stopifnot(nrow(model_matrix) == ncol(Y))
  stopifnot(nrow(beta_mat_init) == nrow(Y))
  stopifnot(ncol(beta_mat_init) == ncol(model_matrix))
  stopifnot(length(dispersions) == nrow(Y))
  stopifnot(dim(offset_matrix) == dim(Y))
  stopifnot(is.null(ridge_penalty) ||
              (is.matrix(ridge_penalty) && ncol(ridge_penalty) == ncol(model_matrix)) ||
              length(ridge_penalty) == ncol(model_matrix))


  if(! is.null(ridge_penalty) && ! is.matrix(ridge_penalty)){
    ridge_target <- attr(ridge_penalty, "target")
    ridge_penalty <- diag(ridge_penalty, nrow = length(ridge_penalty))
    attr(ridge_penalty, "target") <- ridge_target
  }

  apply_ridge <- ! is.null(ridge_penalty)
  n_samples <- ncol(Y)
  if(apply_ridge){
    ridge_penalty_sq <- t(ridge_penalty) %*% ridge_penalty
    ridge_target <- if(is.null(attr(ridge_penalty, "target", TRUE))){
      rep(0, ncol(model_matrix))
    }else{
      attr(ridge_penalty, "target", TRUE)
    }
  }
  result <- list(Beta = matrix(NA, nrow = nrow(Y), ncol = ncol(model_matrix)),
                 iterations = rep(NA, nrow(Y)),
                 deviances = rep(NA, nrow(Y)))

  for(idx in seq_len(nrow(Y))){
    y <- Y[idx, ]
    off <- offset_matrix[idx, ]
    beta_init <- beta_mat_init[idx, ]
    theta <- dispersions[idx]
    if(! apply_ridge){
      res <- optim(par = beta_init, function(beta){
        mu <- exp(model_matrix %*% beta + off)
        compute_gp_deviance_sum(y, mu, theta)
      }, method = "BFGS", control = list(maxit = max_iter))
    }else{
      res <- optim(par = beta_init, function(beta){
        mu <- exp(model_matrix %*% beta + off)
        penalty <- n_samples * t(beta - ridge_target) %*% ridge_penalty_sq %*% (beta - ridge_target)
        compute_gp_deviance_sum(y, mu, theta) + penalty
      }, method = "BFGS", control = list(maxit = max_iter))
    }
    if(res$convergence != 0){
      result$iterations[idx] <- max_iter
      result$Beta[idx, ] <- NA_real_
      result$deviances[idx] <- NA_real_
    }else{
      result$iterations[idx] <- min(res$counts[1], max_iter - 1)
      result$Beta[idx, ] <- res$par
      result$deviances[idx] <- res$value
    }
  }

  result
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
       deviances = matrixStats::rowSums2(Deviance_mat))
}

estimate_betas_group_wise_optimize_helper <- function(y, offset, theta, lower_bound = -30, upper_bound = 30){
  optimize(function(beta){
    sum(dnbinom(y, mu = exp(beta + offset), size = 1/theta, log = TRUE))
  }, lower = lower_bound, upper = upper_bound, maximum = TRUE)$maximum
}

