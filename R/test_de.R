

#' Test for Differential Expression
#'
#' Conduct a quasi-likelihood ratio test for a Gamma-Poisson fit.
#'
#' @inheritParams glm_gp
#' @param fit object of class `glmGamPoi`. Usually the result of calling `glm_gp(data, ...)`
#' @param contrast The contrast to test. Can be a single column name (quoted or as a string)
#'   that is removed from the  full model matrix of `fit`. Or a complex contrast comparing two
#'   or more columns: e.g. `A - B`, `"A - 3 * B"`, `(A + B) / 2 - C` etc. \cr
#'   Only one of `contrast` or `reduced_design` must be specified.
#' @param reduced_design a specification of the reduced design used as a comparison to see what
#'   how much better `fit` describes the data.
#'   Analogous to the `design` parameter in `glm_gp()`, it can be either a `formula`,
#'   a `model.matrix()`, or a `vector`. \cr
#'    Only one of `contrast` or `reduced_design` must be specified.
#' @param pval_adjust_method one of the p-value adjustment method from
#'   [p.adjust.methods]. Default: `"BH"`.
#' @param sort_by the name of the column used to sort the result or `NULL` to return an
#'   unsorted table. Default: `NULL`
#' @param decreasing boolean to decide if the result is sorted increasing or decreasing
#'   order. Default: `FALSE`.
#' @param n_max the maximum number of rows to return. Default: `Inf` which means that all
#'   rows are returned
#'
#' @return a `data.frame` with the following columns
#' \describe{
#'   \item{name}{the rownames of the input data}
#'   \item{pval}{the p-values of the quasi-likelihood ratio test}
#'   \item{adj_pval}{the adjusted p-values returned from [p.adjust()]}
#'   \item{f_statistic}{the F-statistic: \eqn{F = (Dev_full - Dev_red) / (df_1 * disp_ql-shrunken)}}
#'   \item{df1}{the degrees of freedom of the test: `ncol(design) - ncol(reduced_design)`}
#'   \item{df2}{the degrees of freedom of the fit: `ncol(data) - ncol(design) + df_0`}
#'   \item{lfc}{the log2-fold change. If the alternative model is specified by `reduced_design` will
#'    be `NA`.}
#' }
#'
#' @examples
#'   Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
#'   annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
#'                       cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
#'   design <- model.matrix(~ group + cont1 + cont2, data = annot)
#'   fit <- glm_gp(Y, design = design)
#'   res <- test_de(fit, reduced_design = ~ group + cont1, col_data = annot)
#'   head(res)
#'
#'   res2 <- test_de(fit, contrast = cont2)
#'   head(res2)
#'
#' @references
#' * Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012). Detecting differential expression
#'   in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical
#'   Applications in Genetics and Molecular Biology, 11(5).
#'   [https://doi.org/10.1515/1544-6115.1826](https://doi.org/10.1515/1544-6115.1826).
#'
#' @export
test_de <- function(fit,
                    contrast,
                    reduced_design = NULL,
                    pval_adjust_method = "BH", sort_by = NULL, decreasing = FALSE,
                    n_max = Inf,
                    verbose = FALSE){

  if(is.null(reduced_design) == missing(contrast)){
    stop("Please provide either an alternative design (formula or matrix) or a contrast ",
         "(name of a column in fit$model_matrix or a combination of them).")
  }
  if(is.null(fit$overdispersion_shrinkage_list)){
    stop("fit$overdispersion_shrinkage_list is NULL. Run 'glm_gp' with ",
         "'overdispersion_shrinkage = TRUE'.")
  }
  disp_trend <- fit$overdispersion_shrinkage_list$dispersion_trend
  if(! missing(contrast)){
    # e.g. a vector with c(1, 0, -1, 0) for contrast = A - C
    cntrst <- parse_contrast(contrast, levels = colnames(fit$model_matrix),
                             direct_call = FALSE)
    # The modifying matrix of reduced_design has ncol(model_matrix) - 1 columns and rank.
    # The columns are all orthogonal to cntrst.
    # see: https://scicomp.stackexchange.com/a/27835/36204
    # Think about this as a rotation of of the design matrix. The QR decomposition approach
    # has the added benefit that the new columns are all orthogonal to each other, which
    # isn't necessary, but makes fitting more robust
    # The following is a simplified version of edgeR's glmLRT (line 159 in glmfit.R)
    qrc <- qr(cntrst)
    rot <- qr.Q(qrc, complete = TRUE)[,-1,drop=FALSE]
    reduced_design <- fit$model_matrix %*% rot
    lfc <- fit$Beta %*% cntrst / log(2)
  }else{
    lfc <- NA
  }
  if(verbose){message("Fit reduced model")}
  data <- fit$data
  do_on_disk <- is_on_disk.glmGamPoi(fit)
  fit_alt <- glm_gp(data, design = reduced_design,
                    size_factors = fit$size_factors,
                    overdispersion = disp_trend,
                    overdispersion_shrinkage = FALSE,
                    on_disk = do_on_disk)

  if(ncol(fit_alt$model_matrix) >= ncol(fit$model_matrix)){
    stop("The reduced model is as complex (or even more complex) than ",
         "the 'fit' model. The 'reduced_design' should contain fewer terms ",
         "the original 'design'.")
  }

  # Likelihood ratio
  if(verbose){message("Calculate quasi likelihood ratio")}
  lr <- fit_alt$deviances - fit$deviances
  df_test <- ncol(fit$model_matrix) - ncol(fit_alt$model_matrix)
  df_test <- ifelse(df_test == 0, NA, df_test)
  df_fit <- fit$overdispersion_shrinkage_list$ql_df0 + (ncol(data) - ncol(fit_alt$model_matrix))
  f_stat <- lr / df_test / fit$overdispersion_shrinkage_list$ql_disp_shrunken
  pval <- pf(f_stat, df_test, df_fit, lower.tail = FALSE, log.p = FALSE)
  adj_pval <- p.adjust(pval, method = pval_adjust_method)

  names <- rownames(data)
  if(is.null(names)){
    names <- sprintf(paste0("row_%0", floor(log10(nrow(data))), "i"), seq_len(nrow(data)))
  }
  if(verbose){message("Preprare results")}
  res <- data.frame(name = names, pval = pval, adj_pval = adj_pval,
                    f_statistic = f_stat, df1 = df_test, df2 = df_fit,
                    lfc = lfc,
                    stringsAsFactors = FALSE, row.names = NULL)
  res <- if(is.null(sort_by)) {
    res
  }else{
    res[order(res[[sort_by]], decreasing = decreasing), ]
  }
  res[seq_len(min(nrow(res), n_max)), ,drop=FALSE]
}



