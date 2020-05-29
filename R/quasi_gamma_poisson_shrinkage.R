

#' Quasi-Likelihood Ratio test for the Gamma-Poisson
#'
#'
#' @export
gampoi_test_qlr <- function(data, fit,
                            reduced_design = NULL,
                            coef = NULL,
                            col_data = NULL,
                            pval_adjust_method = "BH", sort_by = NULL, decreasing = FALSE,
                            n_max = Inf,
                            verbose = FALSE){

  if(is.null(reduced_design) == is.null(coef)){
    stop("Please provide either an alternative design (formula or matrix) or a coef ",
         "(name of a column in fit$model_matrix).")
  }
  if(is.null(fit$overdispersion_shrinkage_list)){
    stop("fit$overdispersion_shrinkage_list is NULL. Run 'glm_gp' with ",
         "'overdispersion_shrinkage = TRUE'.")
  }
  disp_trend <- fit$overdispersion_shrinkage_list$dispersion_trend
  if(! is.null(coef)){
    reduced_design <- fit$model_matrix[,-coef,drop=FALSE]
  }
  Y <- handle_data_parameter(data, on_disk = NULL)
  if(verbose){message("Fit reduced model")}
  fit_alt <- glm_gp(data, design = reduced_design, col_data = col_data,
                    size_factors = fit$size_factors,
                    overdispersion = disp_trend,
                    overdispersion_shrinkage = FALSE)


  # Likelihood ratio
  if(verbose){message("Calculate quasi likelihood ratio")}
  deviance_full <- DelayedMatrixStats::rowSums2(compute_gp_deviance_residuals_matrix(Y, fit$Mu, disp_trend)^2)
  deviance_alt <- DelayedMatrixStats::rowSums2(compute_gp_deviance_residuals_matrix(Y, fit_alt$Mu, disp_trend)^2)
  lr <- deviance_alt - deviance_full
  df_test <- ncol(fit$model_matrix) - ncol(fit_alt$model_matrix)
  df_fit <- fit$overdispersion_shrinkage_list$ql_df0 + (ncol(Y) - ncol(fit_alt$model_matrix))
  f_stat <- lr / df_test / fit$overdispersion_shrinkage_list$ql_disp_shrunken
  pval <- pf(f_stat, df_test, df_fit, lower.tail = FALSE, log.p = FALSE)
  adj_pval <- p.adjust(pval, method = pval_adjust_method)

  names <- rownames(data)
  if(is.null(names)){
    names <- sprintf(paste0("row_%0", floor(log10(nrow(data))), "i"), seq_len(nrow(data)))
  }
  if(verbose){message("Preprare results")}
  res <- data.frame(name = names, pval = pval, adj_pval = adj_pval, f_statistic = f_stat, df1 = df_test, df2 = df_fit)
  res <- if(is.null(sort_by)) {
    res
  }else{
    res[order(res[[sort_by]], decreasing = decreasing), ]
  }
  res[seq_len(min(nrow(res), n_max)), ,drop=FALSE]
}



shrink_ql_dispersion <- function(disp_est, gene_means,
                                 df, disp_trend = NULL,
                                 npoints = min(length(disp_est), max(0.1 * length(disp_est), 100)), weighted = TRUE,
                                 verbose = FALSE){

  if(verbose){ message("Fit dispersion trend") }
  if(is.null(disp_trend) || isTRUE(disp_trend)){
    disp_trend <- loc_median_fit(gene_means, y = disp_est, npoints = npoints, weighted = weighted)
  }

  if(verbose){ message("Shrink dispersion estimates") }
  # The following equations are used to go between quasi-likelihood and normal representation
  # variance = (gene_means + disp_est * gene_means^2)
  # variance = ql_disp * (gene_means + disp_trend * gene_means^2)
  variances <- gene_means + disp_est * gene_means^2
  ql_disp <- variances / (gene_means + disp_trend * gene_means^2)
  # Lund et al. 2012. Equation 4. Not ideal for very small counts
  # ql_disp <- rowSums2(compute_gp_deviance_residuals_matrix(Y, Mu, disp_trend)^2) / (subsample - ncol(model_matrix))

  var_pr <- variance_prior(ql_disp, df, covariate = gene_means)
  list(dispersion_trend = disp_trend, ql_disp_estimate = ql_disp,
       ql_disp_trend = var_pr$variance0, ql_disp_shrunken = var_pr$var_post,
       ql_df0 = var_pr$df0)
}








#' Estimate the scale and df for a Inverse Chisquare distribution that generate the true gene variances
#'
#' This function implements Smyth's 2004 variance shrinkage. It also supports covariates that are
#' fitted to log(s2) with natural splines. This is based on the 2012 Lund et al. quasi-likelihood
#' paper.
#'
#' @param s2 vector of observed variances. Must not contain \code{0}'s.
#' @param df vector or single number with the degrees of freedom
#' @param covariate a vector with the same length as s2. \code{covariate} is used to regress
#' out the trend in \code{s2}. If \code{covariate = NULL}, it is ignored.
#'
#' @return a list with three elements:
#' \describe{
#'   \item{variance0}{estimate of the scale of the inverse Chisquared distribution. If
#'   covariate is \code{NULL} a single number, otherwise a vector of \code{length(covariate)}}
#'   \item{df0}{estimate of the degrees of freedom of the inverse Chisquared distribution}
#'   \item{var_pos}{the shrunken variance estimates: a combination of \code{s2} and \code{variance0}}
#' }
#'
#' @seealso \code{limma::squeezeVar()}
#'
#' @keywords internal
variance_prior <- function(s2, df, covariate = NULL){
  stopifnot(length(s2) == length(df) || length(df) == 1)
  stopifnot(is.null(covariate) || length(covariate) == length(s2))
  if(any(df <= 0, na.rm=TRUE)){
    stop(paste0("All df must be positive. ", paste0(which(df < 0), collapse=", "), " are not."))
  }
  if(any(s2 <= 0, na.rm=TRUE)){
    stop(paste0("All s2 must be positive. ", paste0(which(s2 < 0), collapse=", "), " are not."))
  }
  if(is.null(covariate)){
    opt_res <- optim(par=c(log_variance0=0, log_df0_inv=0), function(par){
      -sum(df(s2/exp(par[1]), df1=df, df2=exp(par[2]), log=TRUE) - par[1], na.rm=TRUE)
    })

    variance0 <- exp(unname(opt_res$par[1]))
    df0 <- exp(unname(opt_res$par[2]))
  }else{
    # Fit a spline through the log(s2) ~ covariate + 1 to account
    # for remaing trends in s2
    design <- splines::ns(log(covariate), df = 4, intercept = TRUE)
    init_fit <- lm.fit(design, log(s2))
    opt_res <- optim(par=c(betas = init_fit$coefficients, log_df0_inv=0), function(par){
      variance0 <- exp(design %*% par[1:4])
      df0 <- exp(par[5])
      -sum(df(s2/variance0, df1=df, df2=df0, log=TRUE) - log(variance0), na.rm=TRUE)
    }, control = list(maxit = 5000))

    variance0 <- c(exp(design %*% opt_res$par[1:4]))
    df0 <- exp(unname(opt_res$par[5]))
  }

  if(opt_res$convergence != 0){
    warning("Variance prior estimate did not properly converge\n")
    print(opt_res)
  }

  var_post <- (df0 * variance0 + df * s2) / (df0 + df)
  list(variance0 = variance0, df0 = df0, var_post = var_post)
}
