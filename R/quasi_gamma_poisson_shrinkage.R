


#' Shrink the overdispersion estimates
#'
#' Low-level function to shrink a set of overdispersion
#' estimates following the quasi-likelihood and Empirical
#' Bayesian framework.
#'
#' @param disp_est vector of overdispersion estimates
#' @param gene_means vector of average gene expression values. Used
#'   to fit `disp_trend` if that is `NULL`.
#' @param df degrees of freedom for estimating the Empirical Bayesian
#'   variance prior. Can be length 1 or same length as `disp_est` and
#'   `gene_means`.
#' @param disp_trend vector with the dispersion trend. If `NULL` or `TRUE` the
#'   dispersion trend is fitted using a (weighted) local median fit.
#'   Default: `TRUE`.
#' @param ql_disp_trend a logical to indicate if a second abundance trend
#'   using splines is fitted for the quasi-likelihood dispersions.
#'   Default: `NULL` which means that the extra fit is only done if
#'   enough observations are present.
#' @param ... additional parameters for the `loc_median_fit()` function
#' @inheritParams glm_gp
#'
#'
#'
#' @details
#' The function goes through the following steps
#' 1. Fit trend between overdispersion MLE's and the average
#'  gene expression. Per default it uses the `loc_median_fit()`
#'  function.
#' 2. Convert the overdispersion MLE's to quasi-likelihood
#'  dispersion estimates by fixing the trended dispersion as
#'  the "true" dispersion value:
#'  \eqn{disp_ql = (1 + mu * disp_mle) / (1 + mu * disp_trend)}
#' 3. Shrink the quasi-likelihood dispersion estimates using
#'  Empirical Bayesian variance shrinkage (see Smyth 2004).
#'
#' @return the function returns a list with the following elements
#' \describe{
#'   \item{dispersion_trend}{the dispersion trend provided by `disp_trend` or the
#'     local median fit.}
#'   \item{ql_disp_estimate}{the quasi-likelihood dispersion estimates based on
#'     the dispersion trend, `disp_est`, and `gene_means`}
#'   \item{ql_disp_trend}{the `ql_disp_estimate` still might show a trend with
#'     respect to `gene_means`. If `ql_disp_trend = TRUE` a spline is used to
#'     remove this secondary trend. If `ql_disp_trend = TRUE` it corresponds
#'     directly to the dispersion prior}
#'   \item{ql_disp_shrunken}{the shrunken quasi-likelihood dispersion estimates.
#'     They are shrunken towards `ql_disp_trend`.}
#'   \item{ql_df0}{the degrees of freedom of the empirical Bayesian shrinkage.
#'     They correspond to spread of the `ql_disp_estimate`'s}
#' }
#'
#' @examples
#'  Y <- matrix(rnbinom(n = 300 * 4, mu = 6, size = 1/4.2), nrow = 30, ncol = 4)
#'  disps <- sapply(seq_len(nrow(Y)), function(idx){
#'    overdispersion_mle(Y[idx, ])$estimate
#'  })
#'  shrink_list <- overdispersion_shrinkage(disps, rowMeans(Y), df = ncol(Y) - 1,
#'                                          disp_trend = FALSE, ql_disp_trend = FALSE)
#'
#'  plot(rowMeans(Y), shrink_list$ql_disp_estimate)
#'  lines(sort(rowMeans(Y)), shrink_list$ql_disp_trend[order(rowMeans(Y))], col = "red")
#'  points(rowMeans(Y), shrink_list$ql_disp_shrunken, col = "blue", pch = 16, cex = 0.5)
#'
#'
#' @references
#' * Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012). Detecting differential expression
#'   in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical
#'   Applications in Genetics and Molecular Biology, 11(5).
#'   [https://doi.org/10.1515/1544-6115.1826](https://doi.org/10.1515/1544-6115.1826).
#' * Smyth, G. K. (2004). Linear models and empirical bayes methods for assessing differential expression
#'   in microarray experiments. Statistical Applications in Genetics and Molecular Biology, 3(1).
#'   [https://doi.org/10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027)
#'
#' @seealso \code{limma::squeezeVar()}
#' @export
overdispersion_shrinkage <- function(disp_est, gene_means, df,
                                     disp_trend = TRUE,
                                     ql_disp_trend = NULL,
                                     ...,
                                     verbose = FALSE){
  est_value <- ! (is.na(disp_est) | is.na(gene_means) | is.na(df))
  stopifnot(length(disp_est) == length(gene_means))

  if(is.null(disp_trend) || isTRUE(disp_trend)){
    if(verbose){ message("Fit dispersion trend") }
    disp_trend <- rep(NA, length(disp_est))
    disp_trend[est_value] <- loc_median_fit(gene_means[est_value], y = disp_est[est_value], ...)
  }else if(isFALSE(disp_trend)) {
    disp_trend <- rep(NA, length(disp_est))
    disp_trend[est_value] <- rep(mean(disp_est, na.rm = TRUE), sum(est_value))
  }

  if(verbose){ message("Shrink dispersion estimates") }
  # Lund et al. 2012. Equation 4. Not ideal for very small counts
  # ql_disp <- rowSums2(compute_gp_deviance_residuals_matrix(Y, Mu, disp_trend)^2) / df

  # The following equations are used to go between quasi-likelihood and normal representation
  # variance = (gene_means + disp_est * gene_means^2)
  # variance = ql_disp * (gene_means + disp_trend * gene_means^2)
  ql_disp <- (1 + gene_means * disp_est) / (1 + gene_means * disp_trend)

  var_pr <- variance_prior(ql_disp, df, covariate = gene_means, abundance_trend = ql_disp_trend)

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
#'   out the trend in \code{s2}. If \code{covariate = NULL}, it is ignored.
#' @param abundance_trend logical that decides if the additional abundance trend is fit
#'   to the data. If `NULL` the abundance trend is fitted if there are more than 10 observations
#'   and the `covariate` is not `NULL`. Default: `NULL`
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
variance_prior <- function(s2, df, covariate = NULL, abundance_trend = NULL){
  stopifnot(length(s2) == length(df) || length(df) == 1)
  stopifnot(is.null(covariate) || length(covariate) == length(s2))
  if(is.null(covariate) && isTRUE(abundance_trend)) {
    stop("If abundance_trend is true, covariate must not be 'NULL'.")
  }

  if(all(is.na(s2))){
    # This happens if input has zero columns
    return(list(variance0 = rep(NA, length(s2)), df0 = NA, var_post = s2))
  }

  s2_sub <- na.omit(s2)
  if(all(s2_sub == 1)){
    # Happens for example if overdispersion is fixed to 0 -> Poisson GLM
    return(list(variance0 = rep(1, length(s2)), df0 = Inf, var_post = s2))
  }
  if(any(df <= 0, na.rm=TRUE)){
    stop(paste0("All df must be positive. df ", paste0(which(df <= 0), collapse=", "), " are not."))
  }
  if(any(s2 <= 0, na.rm=TRUE)){
    stop(paste0("All s2 must be positive. s2 ", paste0(which(s2 <= 0), collapse=", "), " are not."))
  }

  # Fit an intercept to s2
  opt_res <- optim(par=c(log_variance0=0, log_df0_inv=0), function(par){
    -sum(df(s2_sub/exp(par[1]), df1=df, df2=exp(par[2]), log=TRUE) - par[1], na.rm=TRUE)
  })
  variance0 <- rep(exp(unname(opt_res$par[1])), times = length(s2_sub))
  df0 <- exp(unname(opt_res$par[2]))

  if(is.numeric(abundance_trend)){
    stop("Numeric abundance trend is not supported.")
  }else if(is.null(covariate) || is.null(abundance_trend) || ! abundance_trend || length(s2_sub) < 10 || sum(covariate > 1e-8, na.rm = TRUE) < 10) {
    if(length(s2_sub) < 10 && isTRUE(abundance_trend)){
      warning("abundance_trend was set to 'TRUE', however there were not enough observations (length(s2) = ",
              length(s2_sub), " < 10) to fit the trend.")
    }
    # Do nothing else to data
  }else{
    # Fit a trend!
    # Fit a spline through the log(s2_sub) ~ covariate + 1 to account
    # for remaing trends in s2_sub
    covariate_sub <- covariate[! is.na(s2)]
    log_covariate <- log(covariate_sub)
    not_zero <- covariate_sub > 1e-8
    ra <- range(log_covariate[not_zero])
    knots <- ra[1] + c(1/3, 2/3) * (ra[2] - ra[1])
    tryCatch({
      design <- splines::ns(log_covariate[not_zero], df = 4, knots = knots, intercept = TRUE)
      init_fit <- lm.fit(design, log(s2_sub[not_zero]))
      opt_res <- optim(par=c(betas = init_fit$coefficients, log_df0_inv=0), function(par){
        variance0 <- exp(design %*% par[1:4])
        df0 <- exp(par[5])
        -sum(df(s2_sub[not_zero]/variance0, df1=df, df2=df0, log=TRUE) - log(variance0), na.rm=TRUE)
      }, control = list(maxit = 5000))
      variance0[not_zero] <- c(exp(design %*% opt_res$par[1:4]))
      df0 <- exp(unname(opt_res$par[5]))
    }, error = function(e){
      warning("Problem fitting the abundance trend to the quasi-likelihood dispersion. ",
              "Falling back to non-trended variance0 estimation.")
    })
  }

  if(opt_res$convergence != 0){
    warning("Variance prior estimate did not properly converge\n")
  }

  var_post <- rep(NA, length(s2))
  var_post[! is.na(s2)] <- (df0 * variance0 + df * s2_sub) / (df0 + df)
  list(variance0 = variance0, df0 = df0, var_post = var_post)
}
