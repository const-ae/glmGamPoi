

#' Test for Differential Expression
#'
#' Conduct a quasi-likelihood ratio test for a Gamma-Poisson fit.
#'
#' @inheritParams glm_gp
#' @param fit object of class `glmGamPoi`. Usually the result of calling `glm_gp(data, ...)`
#' @param contrast The contrast to test. Can be a single column name (quoted or as a string)
#'   that is removed from the  full model matrix of `fit`. Or a complex contrast comparing two
#'   or more columns: e.g. `A - B`, `"A - 3 * B"`, `(A + B) / 2 - C` etc. For complicated
#'   experimental design that involved nested conditions, you specify the condition level to compare
#'   using the `fact()` helper function. \cr
#'   Only one of `contrast` or `reduced_design` must be specified.
#' @param reduced_design a specification of the reduced design used as a comparison to see what
#'   how much better `fit` describes the data.
#'   Analogous to the `design` parameter in `glm_gp()`, it can be either a `formula`,
#'   a `model.matrix()`, or a `vector`. \cr
#'    Only one of `contrast` or `reduced_design` must be specified.
#' @param full_design option to specify an alternative `full_design` that can differ from
#'   `fit$model_matrix`. Can be a `formula` or a `matrix`. Default: `fit$model_matrix`
#' @param subset_to a vector with the same length as `ncol(fit$data)` or  an expression
#'   that evaluates to such a vector. The expression can reference columns from `colData(fit$data)`.
#'   A typical use case in single cell analysis would be to subset to a specific cell type
#'   (e.g. `subset_to = cell_type == "T-cells"`). Note that if this argument is set a new
#'   the model for the `full_design` is re-fit.\cr
#'   Default: `NULL` which means that the data is not subset.
#' @param pseudobulk_by *DEPRECTATED*, please use the [pseudobulk] function instead. \cr
#'   A vector with the same length as `ncol(fit$data)` that is used to
#'   split the columns into different groups (calls [split()]). `pseudobulk_by` can also be an
#'   expression that evaluates to a vector. The expression can reference columns from `colData(fit$data)`. \cr
#'   The counts are summed across the groups
#'   to create "pseudobulk" samples. This is typically used in single cell analysis if the cells come
#'   from different samples to get a proper estimate of the differences. This is particularly powerful
#'   in combination with the `subset_to` parameter to analyze differences between samples for
#'   subgroups of cells. Note that this does a fresh fit for both the full and the reduced design.
#'   Default: `NULL` which means that the data is not aggregated.
#' @param pval_adjust_method one of the p-value adjustment method from
#'   [p.adjust.methods]. Default: `"BH"`.
#' @param sort_by the name of the column or an expression used to sort the result. If `sort_by` is `NULL`
#'   the table is not sorted. Default: `NULL`
#' @param decreasing boolean to decide if the result is sorted increasing or decreasing
#'   order. Default: `FALSE`.
#' @param n_max the maximum number of rows to return. Default: `Inf` which means that all
#'   rows are returned
#' @param max_lfc set the maximum absolute log fold change that is returned. Large log fold changes
#'   occur for lowly expressed genes because the ratio of two small numbers can be impractically large. For example, limiting
#'   the range of log fold changes can clarify the patterns in a volcano plot. Default: `10` which
#'   corresponds to a thousand-fold (2^10) increase in expression.
#'
#' @details
#' The `fact()` helper function simplifies the specification of a contrast for complex experimental designs.
#' Instead of working out which combination of coefficients corresponds to a research question,
#' you can simply specify the two conditions that you want to compare.
#'
#' You can only call the `fact` function inside the `contrast` argument. The arguments are the selected factor levels
#' for each covariate. To compare two conditions, simply subtract the two `fact` calls. Internally, the package
#' calls [model.matrix] using the provided values and the original formula from the fit to produce a vector.
#' Subtracting two of these vectors produces a contrast vector. Missing covariates are filled with the first factor level
#' or zero for numerical covariates.
#'
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
#'  # Make Data
#'  Y <- matrix(rnbinom(n = 30 * 100, mu = 4, size = 0.3), nrow = 30, ncol = 100)
#'  annot <- data.frame(mouse = sample(LETTERS[1:6], size = 100, replace = TRUE),
#'                      celltype = sample(c("Tcell", "Bcell", "Macrophages"), size = 100, replace = TRUE),
#'                      cont1 = rnorm(100), cont2 = rnorm(100, mean = 30))
#'  annot$condition <- ifelse(annot$mouse %in% c("A", "B", "C"), "ctrl", "treated")
#'  head(annot)
#'  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = annot)
#'
#'  # Fit model
#'  fit <- glm_gp(se, design = ~ condition + celltype + cont1 + cont2)
#'  # Test with reduced design
#'  res <- test_de(fit, reduced_design = ~ celltype + cont1 + cont2)
#'  head(res)
#'  # Test with contrast argument, the results are identical
#'  res2 <- test_de(fit, contrast = conditiontreated)
#'  head(res2)
#'  # Test with explicit specification of the conditions
#'  # The results are still identical
#'  res3 <- test_de(fit, contrast = fact(condition = "treated", celltype = "Bcell") -
#'                                     fact(condition = "ctrl", celltype = "Bcell"))
#'  head(res3)
#'
#'
#'  # The column names of fit$Beta are valid variables in the contrast argument
#'  colnames(fit$Beta)
#'  # You can also have more complex contrasts:
#'  # the following compares cont1 vs cont2:
#'  test_de(fit, cont1 - cont2, n_max = 4)
#'  # You can also sort the output
#'  test_de(fit, cont1 - cont2, n_max = 4,
#'          sort_by = "pval")
#'  test_de(fit, cont1 - cont2, n_max = 4,
#'          sort_by = - abs(f_statistic))
#'
#'  # If the data has multiple samples, it is a good
#'  # idea to aggregate the cell counts by samples.
#'  # This is called forming a "pseudobulk".
#'  se_reduced <- pseudobulk(se, group_by = vars(mouse, condition, celltype),
#'                           cont1 = mean(cont1), cont2 = min(cont2))
#'  fit_reduced <- glm_gp(se_reduced, design = ~ condition + celltype)
#'  test_de(fit_reduced, contrast = "conditiontreated", n_max = 4)
#'  test_de(fit_reduced, contrast = fact(condition = "treated", celltype = "Macrophages") - fact(condition = "ctrl", celltype = "Macrophages"),
#'          n_max = 4)
#'
#'
#' @references
#' * Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012). Detecting differential expression
#'   in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical
#'   Applications in Genetics and Molecular Biology, 11(5).
#'   [https://doi.org/10.1515/1544-6115.1826](https://doi.org/10.1515/1544-6115.1826).
#'
#' @seealso [glm_gp()]
#'
#' @export
test_de <- function(fit,
                    contrast,
                    reduced_design = NULL,
                    full_design = fit$model_matrix,
                    subset_to = NULL, pseudobulk_by = NULL,
                    pval_adjust_method = "BH", sort_by = NULL,
                    decreasing = FALSE, n_max = Inf,
                    max_lfc = 10,
                    verbose = FALSE){
  # Capture all NSE variables
  subset_to_capture <- substitute(subset_to)
  pseudobulk_by_capture <- substitute(pseudobulk_by)
  sort_by_capture <- substitute(sort_by)
  test_de_q(fit, contrast = {{contrast}}, reduced_design = reduced_design, full_design = full_design,
            subset_to = subset_to_capture, pseudobulk_by = pseudobulk_by_capture,
            pval_adjust_method = pval_adjust_method, sort_by = sort_by_capture,
            decreasing = decreasing, n_max = n_max, max_lfc = max_lfc,
            verbose = verbose,
            env = parent.frame())
}

# This function is necessary to handle the NSE stuff
# for an explanation see:
# http://adv-r.had.co.nz/Computing-on-the-language.html
test_de_q <- function(fit,
                      contrast,
                      reduced_design = NULL,
                      full_design = fit$model_matrix,
                      subset_to = NULL, pseudobulk_by = NULL,
                      pval_adjust_method = "BH", sort_by = NULL,
                      decreasing = FALSE, n_max = Inf, max_lfc = 10,
                      verbose = FALSE,
                      env = parent.frame()){

  ridge_penalty <- fit$ridge_penalty
  if(! is.matrix(full_design) || length(full_design) != length(fit$model_matrix) || ! all(full_design == fit$model_matrix) ){
    full_design <- handle_design_parameter(full_design, fit$data, SummarizedExperiment::colData(fit$data), NULL)$model_matrix
    ridge_penalty <- NULL
  }
  if(is.null(reduced_design) == rlang::quo_is_missing(rlang::enquo(contrast))){
    stop("Please provide either an alternative design (formula or matrix) or a contrast ",
         "(name of a column in fit$model_matrix or a combination of them).")
  }

  pseudobulk_by_e <- eval_with_q(pseudobulk_by, SummarizedExperiment::colData(fit$data), env = env)
  subset_to_e <- eval_with_q(subset_to, SummarizedExperiment::colData(fit$data), env = env)
  if(! is.null(pseudobulk_by_e)){
    if(verbose) { message("Aggregate counts of columns that belong to the same sample") }
    # This method aggregates the data
    # then does a fresh fit for the full model
    # then calls test_de() with the reduced dataset
    return(test_pseudobulk_q(fit$data, design = full_design,
                             aggregate_cells_by = pseudobulk_by,
                             contrast = {{contrast}},
                             reduced_design = reduced_design,
                             ridge_penalty = ridge_penalty,
                             subset_to = subset_to_e,
                             pval_adjust_method = pval_adjust_method, sort_by = sort_by,
                             decreasing = decreasing, n_max = n_max,
                             verbose = verbose, env = env))
  }

  if(! is.null(subset_to_e)){
    if(verbose) { message("Subset the dataset, thus a new model is fitted!") }
    # Only run test on a subset of things:
    # Redo fit, but keep dispersion
    data_subset <- fit$data[,subset_to_e,drop=FALSE]

    model_matrix_subset <- full_design[subset_to_e, ,drop=FALSE]
    size_factor_subset <- fit$size_factors[subset_to_e]

    fit_subset <- glm_gp(data_subset, design = model_matrix_subset, size_factors = size_factor_subset,
                  overdispersion = fit$overdispersions, on_disk = is_on_disk.glmGamPoi(fit), verbose = verbose)
    test_res <- test_de_q(fit_subset, contrast = {{contrast}}, reduced_design = reduced_design,
                          subset_to = NULL, pseudobulk_by = NULL,
                          pval_adjust_method = pval_adjust_method, sort_by = sort_by,
                          decreasing = decreasing, n_max = n_max, max_lfc = max_lfc,
                          verbose = verbose, env = env)
    return(test_res)
  }
  if(is.null(fit$overdispersion_shrinkage_list)){
    stop("fit$overdispersion_shrinkage_list is NULL. Run 'glm_gp' with ",
         "'overdispersion_shrinkage = TRUE'.")
  }
  disp_trend <- fit$overdispersion_shrinkage_list$dispersion_trend
  if(!rlang::quo_is_missing(rlang::enquo(contrast))){
    # e.g. a vector with c(1, 0, -1, 0) for contrast = A - C
    cntrst <- parse_contrast({{contrast}}, coefficient_names = colnames(fit$model_matrix), formula = fit$design_formula)
    cntrst <- as.matrix(cntrst)
    if(nrow(cntrst) != ncol(fit$model_matrix)){
      stop("The length of the contrast vector does not match the number of coefficients in the model (",
           ncol(fit$model_matrix), ")\n", format_matrix(cntrst))
    }
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
    reduced_ridge_penalty <- if(is.null(ridge_penalty)){
      NULL
    }else if(is.matrix(ridge_penalty)){
      ridge_penalty %*% rot
    }else{
      diag(ridge_penalty, nrow = length(ridge_penalty)) %*% rot
    }
    lfc <- fit$Beta %*% cntrst / log(2)
    lfc[lfc < -max_lfc] <- -max_lfc
    lfc[lfc > max_lfc] <- max_lfc
  }else{
    # Convert the formula to model matrix
    reduced_design <- handle_design_parameter(reduced_design,  fit$data,
                                   SummarizedExperiment::colData(fit$data), reference_level = NULL)$model_matrix
    if(ncol(reduced_design) >= ncol(full_design)){
      stop("The reduced model is as complex (or even more complex) than ",
           "the 'fit' model. The 'reduced_design' should contain fewer terms ",
           "the original 'design'.")
    }

    # Check if full_model is nested in full_design
    rot <- solve_lm_for_B(reduced_design, full_design)
    if(any(abs(reduced_design - full_design %*% rot) > 1e-10)){
      warning("Although the 'reduced_design' matrix has fewer columns than ",
              "'fit$model_matrix', it appears that the 'reduced_design' is not ",
              "nested in the 'fit$model_matrix'. Accordingly, the results of the ",
              "statistical test will be unreliable.")
    }
    reduced_ridge_penalty <- if(is.null(ridge_penalty)){
      NULL
    }else if(is.matrix(ridge_penalty)){
      ridge_penalty %*% rot
    }else{
      diag(ridge_penalty, nrow = length(ridge_penalty)) %*% rot
    }
    lfc <- NA
  }
  if(verbose){message("Fit reduced model")}
  data <- fit$data
  do_on_disk <- is_on_disk.glmGamPoi(fit)
  fit_alt <- glm_gp(data, design = reduced_design,
                    size_factors = FALSE, # size factors are already in offset
                    offset = fit$Offset,
                    overdispersion = disp_trend,
                    overdispersion_shrinkage = FALSE,
                    ridge_penalty = reduced_ridge_penalty,
                    on_disk = do_on_disk)



  # Likelihood ratio
  if(verbose){message("Calculate quasi likelihood ratio")}
  lr <- (fit_alt$deviances - fit$deviances)
  df_test <- ncol(fit$model_matrix) - ncol(fit_alt$model_matrix)
  df_test <- ifelse(df_test == 0, NA, df_test)
  df_fit <- fit$overdispersion_shrinkage_list$ql_df0 + (ncol(data) - ncol(fit$model_matrix))
  f_stat <- lr / df_test / fit$overdispersion_shrinkage_list$ql_disp_shrunken
  pval <- pf(f_stat, df_test, df_fit, lower.tail = FALSE, log.p = FALSE)
  adj_pval <- p.adjust(pval, method = pval_adjust_method)

  names <- rownames(data)
  if(is.null(names)){
    names <- sprintf(paste0("row_%0", floor(log10(nrow(data))), "i"), seq_len(nrow(data)))
  }
  if(verbose){message("Prepare results")}
  res <- data.frame(name = names, pval = pval, adj_pval = adj_pval,
                    f_statistic = f_stat, df1 = df_test, df2 = df_fit,
                    lfc = lfc,
                    stringsAsFactors = FALSE, row.names = NULL)
  sort_by_e <- eval_with_q(sort_by, res, env = env)
  res <- if(is.null(sort_by_e)) {
    res
  }else{
    res[order(sort_by_e, decreasing = decreasing), ]
  }
  res[seq_len(min(nrow(res), n_max)), ,drop=FALSE]
}



