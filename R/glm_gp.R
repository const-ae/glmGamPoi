

#' Fit a Gamma-Poisson Generalized Linear Model
#'
#' This function provides a simple to use interface to fit Gamma-Poisson generalized
#' linear models. It works equally well for small scale (a single model) and large scale data
#' (e.g. thousands of rows and columns, potentially stored on disk). The function
#' automatically determines the appropriate size factors for each sample and efficiently
#' finds the best overdispersion parameter for each gene.
#'
#' @param data any matrix-like object (e.g. [matrix], [DelayedArray], [HDF5Matrix]) or
#'   anything that can be cast to a [SummarizedExperiment] (e.g. `MSnSet`, `eSet` etc.) with
#'   one column per sample and row per gene.
#' @param design a specification of the experimental design used to fit the Gamma-Poisson GLM.
#'   It can be a [model.matrix()] with one row for each sample and one column for each
#'   coefficient. \cr
#'   Alternatively, `design` can be a `formula`. The entries in the
#'   formula can refer to global objects, columns in the `col_data` parameter, or the `colData(data)`
#'   of `data` if it is a `SummarizedExperiment`. \cr
#'   The third option is that `design` is a vector where each element specifies to which
#'   condition a sample belongs. \cr
#'   Default: `design = ~ 1`, which means that all samples are treated as if they belong to the
#'   same condition. Note that this is the fasted option.
#' @param col_data a dataframe with one row for each sample in `data`. Default: `NULL`.
#' @param reference_level a single string that specifies which level is used as reference
#'   when the model matrix is created. The reference level becomes the intercept and all
#'   other coefficients are calculated with respect to the `reference_level`.
#'   Default: `NULL`.
#' @param offset Constant offset in the model in addition to `log(size_factors)`. It can
#'   either be a single number, a vector of length `ncol(data)` or a matrix with the
#'   same dimensions as `dim(data)`. Note that if data is a [DelayedArray] or [HDF5Matrix],
#'   `offset` must be as well. Default: `0`.
#' @param size_factors in large scale experiments, each sample is typically of different size
#'   (for example different sequencing depths). A size factor is an internal mechanism of GLMs to
#'   correct for this effect.\cr
#'   `size_factors` is either a numeric vector with positive entries that has the same lengths as columns in the data
#'   that specifies the size factors that are used.
#'   Or it can be a string that species the method that is used to estimate the size factors
#'   (one of \code{c("normed_sum", "deconvolution", "poscounts")}).
#'   Note that \code{"normed_sum"} and \code{"poscounts"} are fairly
#'   simple methods and can lead to suboptimal results. For the best performance, I recommend to use
#'   `size_factors = "deconvolution"` which calls `scran::calculateSumFactors()`. However, you need
#'   to separately install the `scran` package from Bioconductor for this method to work.
#'   Also note that `size_factors = 1` and `size_factors = FALSE` are equivalent. If only a single gene is given,
#'   no size factor is estimated (ie. `size_factors = 1`). Default: `"normed_sum"`.
#' @param overdispersion the simplest count model is the Poisson model. However, the Poisson model
#'   assumes that \eqn{variance = mean}. For many applications this is too rigid and the Gamma-Poisson
#'   allows a more flexible mean-variance relation (\eqn{variance = mean + mean^2 * overdispersion}). \cr
#'   `overdispersion` can either be
#'   \itemize{
#'      \item a single boolean that indicates if an overdispersion is estimated for each gene.
#'      \item a numeric vector of length `nrow(data)` fixing the overdispersion to those values.
#'      \item the string `"global"` to indicate that one dispersion is fit across all genes.
#'   }
#'   Note that `overdispersion = 0` and `overdispersion = FALSE` are equivalent and both reduce
#'   the Gamma-Poisson to the classical Poisson model. Default: `TRUE`.
#' @param overdispersion_shrinkage the overdispersion can be difficult to estimate with few replicates. To
#'   improve the overdispersion estimates, we can share information across genes and shrink each individual
#'   overdispersion estimate towards a global overdispersion estimate. Empirical studies show however that
#'   the overdispersion varies based on the mean expression level (lower expression level => higher
#'   dispersion). If `overdispersion_shrinkage = TRUE`, a median trend of dispersion and expression level is
#'   fit and used to estimate the variances of a quasi Gamma Poisson model (Lund et al. 2012). Default: `TRUE`.
#' @param do_cox_reid_adjustment the classical maximum likelihood estimator of the `overdisperion` is biased
#'   towards small values. McCarthy _et al._ (2012) showed that it is preferable to optimize the Cox-Reid
#'   adjusted profile likelihood.\cr
#'   `do_cox_reid_adjustment` can be either be `TRUE` or `FALSE` to indicate if the adjustment is
#'   added during the optimization of the `overdispersion` parameter. Default: `TRUE`.
#' @param subsample the estimation of the overdispersion is the slowest step when fitting
#'   a Gamma-Poisson GLM. For datasets with many samples, the estimation can be considerably sped up
#'   without loosing much precision by fitting the overdispersion only on a random subset of the samples.
#'   Default: `FALSE` which means that the data is not subsampled. If set to `TRUE`, at most 1,000 samples
#'   are considered. Otherwise the parameter just specifies the number of samples that are considered
#'   for each gene to estimate the overdispersion.
#' @param on_disk a boolean that indicates if the dataset is loaded into memory or if it is kept on disk
#'   to reduce the memory usage. Processing in memory can be significantly faster than on disk.
#'   Default: `NULL` which means that the data is only processed in memory if `data` is an in-memory
#'   data structure.
#' @param verbose a boolean that indicates if information about the individual steps are printed
#'   while fitting the GLM. Default: `FALSE`.
#'
#'
#' @details
#' The method follows the following steps:
#'
#' 1. The size factors are estimated.\cr
#'    If `size_factors = "normed_sum"`, the column-sum for each cell is calculated and the resulting
#'    size factors are normalized so that their geometric mean is `1`. If `size_factors = "poscounts"`,
#'    a slightly adapted version of the procedure proposed by Anders and Huber (2010) in
#'    equation (5) is used. To handle the large number of zeros the geometric means are calculated for
#'    \eqn{Y + 0.5} and ignored during the calculation of the median. Columns with all zeros get a
#'    default size factor of \eqn{0.001}. If `size_factors = "deconvolution"`, the method
#'    `scran::calculateSumFactors()` is called.
#' 2. The dispersion estimates are initialized based on the moments of each row of \eqn{Y}.
#' 3. The coefficients of the model are estimated.\cr
#'    If all samples belong to the same condition (i.e. `design = ~ 1`), the betas are estimated using
#'    a quick Newton-Raphson algorithm. This is similar to the behavior of `edgeR`. For more complex
#'    designs, the general Fisher-scoring algorithm is used. Here, the code is based on a fork  of the
#'    internal function `fitBeta()` from `DESeq2`. It does however contain some modification to make
#'    the fit more robust and faster.
#' 4. The mean for each gene and sample is calculate.\cr
#'    Note that this step can be very IO intensive if `data` is or contains a DelayedArray.
#' 5. The overdispersion is estimated.\cr
#'    The classical method for estimating the overdispersion for each gene is to maximize the
#'    Gamma-Poisson log-likelihood by iterating over each count and summing the the corresponding
#'    log-likelihood. It is however, much more efficient
#'    for genes with many small counts to work on the contingency table of the counts. Originally, this
#'    approach had already been used by Anscombe (1950). In this package, I have implemented an
#'    extension of their method that can handle general offsets.\cr
#'    See also [overdispersion_mle()].
#'  6. The beta coefficients are estimated once more with the updated overdispersion estimates
#'  7. The mean for each gene and sample is calculated again.
#'
#' This method can handle not just in memory data, but also data stored on disk. This is essential for
#' large scale datasets with thousands of samples, as they sometimes encountered in modern single-cell
#' RNA-seq analysis. `glmGamPoi` relies on the `DelayedArray` and `beachmat` package to efficiently
#' implement the access to the on-disk data.
#'
#' @return
#' The method returns a list with the following elements:
#' \describe{
#'   \item{`Beta`}{a matrix with dimensions `nrow(data) x n_coefficients` where `n_coefficients` is
#'   based on the `design` argument. It contains the estimated coefficients for each gene.}
#'   \item{`overdispersions`}{a vector with length `nrow(data)`. The overdispersion parameter for each
#'   gene. It describes how much more the counts vary than one would expect according to the Poisson
#'   model.}
#'   \item{`overdispersion_shrinkage_list`}{a list with additional information from the quasi-likelihood
#'   shrinkage. For details see [overdispersion_shrinkage()].}
#'   \item{`deviances`}{a vector with the deviance of the fit for each row. The deviance is a measure
#'   how well the data is fit by the model. A smaller deviance means a better fit.}
#'   \item{`Mu`}{a matrix with the same dimensions as `dim(data)`. If the calculation happened on
#'   disk, than `Mu` is a `HDF5Matrix`. It contains the estimated mean value for each gene and
#'   sample.}
#'   \item{`size_factors`}{a vector with length `ncol(data)`. The size factors are the inferred
#'   correction factors for different sizes of each sample. They are also sometimes called the
#'   exposure factor.}
#'   \item{`Offset`}{a matrix with the same dimensions as `dim(data)`. If the calculation happened on
#'   disk, than `Offset` is a `HDF5Matrix`. It contains the `log(size_factors) + offset` from the
#'   function call.}
#'   \item{`data`}{a `SummarizedExperiment` that contains the input counts and the `col_data`}
#'   \item{`model_matrix`}{a matrix with dimensions `ncol(data) x n_coefficients`. It is build based
#'   on the `design` argument.}
#' }
#'
#' @examples
#'  set.seed(1)
#'  # The simplest example
#'  y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
#'  c(glm_gp(y, size_factors = FALSE))
#'
#'  # Fitting a whole matrix
#'  model_matrix <- cbind(1, rnorm(5))
#'  true_Beta <- cbind(rnorm(n = 30), rnorm(n = 30, mean = 3))
#'  sf <- exp(rnorm(n = 5, mean = 0.7))
#'  model_matrix
#'  Y <- matrix(rnbinom(n = 30 * 5, mu = sf * exp(true_Beta %*% t(model_matrix)), size = 1/2.4),
#'              nrow = 30, ncol = 5)
#'
#'  fit <- glm_gp(Y, design = model_matrix, size_factors = sf, verbose = TRUE)
#'  summary(fit)
#'
#'  # Fitting a model with covariates
#'  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
#'  city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
#'  age = rnorm(n = 50, mean = 40, sd = 15))
#'  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
#'  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
#'  summary(fit)
#'
#'
#'
#' @seealso [overdispersion_mle()] and [overdispersion_shrinkage()] for the internal functions that do the
#'   work. For differential expression analysis, see [test_de()].
#'
#' @references
#'   * McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression analysis of multifactor
#'   RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), 4288–4297.
#'   [https://doi.org/10.1093/nar/gks042](https://doi.org/10.1093/nar/gks042).
#'   * Anders Simon, & Huber Wolfgang. (2010). Differential expression analysis for sequence count data.
#'   Genome Biology. [https://doi.org/10.1016/j.jcf.2018.05.006](https://doi.org/10.1016/j.jcf.2018.05.006).
#'   * Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and  dispersion
#'   for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.
#'   [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
#'   * Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2009). edgeR: A Bioconductor package for differential
#'   expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140.
#'   [https://doi.org/10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).
#'   * Lun ATL, Pagès H, Smith ML (2018). “beachmat: A Bioconductor C++ API for accessing high-throughput
#'   biological data from a variety of R matrix types.” PLoS Comput. Biol., 14(5), e1006135. doi:
#'   [10.1371/journal.pcbi.1006135.](https://doi.org/10.1371/journal.pcbi.1006135).
#'   * Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012). Detecting differential expression
#'   in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. Statistical
#'   Applications in Genetics and Molecular Biology, 11(5).
#'   [https://doi.org/10.1515/1544-6115.1826](https://doi.org/10.1515/1544-6115.1826).
#'   * Lun ATL, Bach K and Marioni JC (2016). Pooling across cells to normalize single-cell RNA sequencing
#'   data with many zero counts. Genome Biol. 17:75
#'   [https://doi.org/10.1186/s13059-016-0947-7](https://doi.org/10.1186/s13059-016-0947-7)
#'
#' @export
glm_gp <- function(data,
                   design = ~ 1,
                   col_data = NULL,
                   reference_level = NULL,
                   offset = 0,
                   size_factors = c("normed_sum", "deconvolution", "poscounts"),
                   overdispersion = TRUE,
                   overdispersion_shrinkage = TRUE,
                   ridge_penalty = 0,
                   do_cox_reid_adjustment = TRUE,
                   subsample = FALSE,
                   on_disk = NULL,
                   verbose = FALSE){

  # Validate `data`
  if(inherits(data, "formula")){
    if(length(design) != 2 || design != ~ 1){
      stop("If the first argument is already a formula, the second argument must not be set. Please call this function like this:\n",
           "'glm_gp(data = mat, design = ~ a + b + c, ...)'", call. = FALSE)
    }
    extr <- extract_data_from_formula(data, col_data, parent.frame())
    data <- extr$data
    design <- extr$design
  }
  if(is.vector(data)){
    data <- matrix(data, nrow = 1)
  }
  data_mat <- handle_data_parameter(data, on_disk)

  # Convert the formula to a model_matrix
  col_data <- get_col_data(data, col_data)
  des <- handle_design_parameter(design, data, col_data, reference_level)

  # Call glm_gp_impl()
  res <- glm_gp_impl(data_mat,
              model_matrix = des$model_matrix,
              offset = offset,
              size_factors = size_factors,
              overdispersion = overdispersion,
              overdispersion_shrinkage = overdispersion_shrinkage,
              ridge_penalty = ridge_penalty,
              do_cox_reid_adjustment = do_cox_reid_adjustment,
              subsample = subsample,
              verbose = verbose)
  # Make sure that the output is nice and beautiful
  rownames(data_mat) <- rownames(data)
  colnames(data_mat) <- colnames(data)
  res$data <- SummarizedExperiment::SummarizedExperiment(list(counts = data_mat),
                                                         colData = col_data)
  res$model_matrix <- des$model_matrix
  if(is.null(colnames(res$model_matrix))){
    colnames(res$model_matrix) <- paste0("Coef_", seq_len(ncol(res$model_matrix)))
  }
  res$design_formula <- des$design_formula
  colnames(res$Beta) <- colnames(res$model_matrix)
  rownames(res$Beta) <- rownames(data)
  if(! is.null(res$ridge_penalty)){
    names(res$ridge_penalty) <- colnames(res$model_matrix)
  }
  rownames(res$Mu) <- rownames(data)
  colnames(res$Mu) <- colnames(data)
  rownames(res$Offset) <- rownames(data)
  colnames(res$Offset) <- colnames(data)
  names(res$overdispersions) <- rownames(data)
  names(res$deviances) <- rownames(data)
  names(res$size_factors) <- colnames(data)

  class(res) <- "glmGamPoi"
  res
}


handle_data_parameter <- function(data, on_disk){
  if(is.matrix(data)){
    if(! is.numeric(data)){
      stop("The data argument must consist of numeric values and not of ", mode(data), " values")
    }
    if(is.null(on_disk) || isFALSE(on_disk)){
      data_mat <- data
    }else if(isTRUE(on_disk)){
      data_mat <- HDF5Array::writeHDF5Array(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "DelayedArray")){
    if(is.null(on_disk) || isTRUE(on_disk)){
      data_mat <- data
    }else if(isFALSE(on_disk)){
      data_mat <- as.matrix(data)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(data, "SummarizedExperiment")){
    data_mat <- handle_data_parameter(SummarizedExperiment::assay(data), on_disk)
  }else if(canCoerce(data, "SummarizedExperiment")){
    se <- as(data, "SummarizedExperiment")
    data_mat <- handle_data_parameter(SummarizedExperiment::assay(se), on_disk)
  }else if(is(data, "dgCMatrix") || is(data, "dgTMatrix")) {
    stop("glmGamPoi does not yet support sparse input data of type '", class(data),"'. ",
         "If your data fits into memory, please cast it to a dense matrix with ",
         "'as.matrix(data)' or if it is too big, convert it to an on-disk dataset ",
         "with the 'HDF5Array' package. ")
  }else{
    stop("Cannot handle data of class '", class(data), "'.",
         "It must be of type dense matrix object (i.e., a base matrix or DelayedArray),",
         " or a container for such a matrix (for example: SummarizedExperiment).")
  }
  data_mat
}


get_col_data <- function(data, col_data){
  if(is.null(col_data)){
    col_data <- as.data.frame(matrix(numeric(0), nrow=ncol(data)))
  }
  if(is(data, "SummarizedExperiment")){
    cbind(col_data, SummarizedExperiment::colData(data))
  }else{
    col_data
  }
}


handle_design_parameter <- function(design, data, col_data, reference_level){
  n_samples <- ncol(data)

  # Handle the design parameter
  if(is.matrix(design)){
    if(! is.null(reference_level)){
      stop("Cannot specify `design` as a matrix and `reference_level` != NULL.\n",
           "Either provide `design` as a formula and specify `reference_level` or ",
           "make sure that the correct reference level is used when the matrix is created ",
           "for example with `stats::relevel()`.")
    }
    model_matrix <- design
    design_formula <- NULL
  }else if((is.vector(design) || is.factor(design))){
    if(length(design) != n_samples){
      if(length(design) == 1 && design == 1){
        stop("The specified design vector length (", length(design), ") does not match ",
             "the number of samples: ", n_samples, "\n",
             "Did you maybe mean: `design = ~ 1`?")
      }else{
        stop("The specified design vector length (", length(design), ") does not match ",
             "the number of samples: ", n_samples)
      }
    }
    tmp <- convert_chr_vec_to_model_matrix(design, reference_level)
    model_matrix <- tmp$model_matrix
    design_formula <- tmp$formula
    attr(design_formula, "constructed_from") <- "vector"
  }else if(inherits(design,"formula")){
    col_data <- relevel_columns_to_reference_level(col_data, reference_level)
    tmp <- convert_formula_to_model_matrix(design, col_data)
    model_matrix <- tmp$model_matrix
    design_formula <- tmp$formula
    attr(design_formula, "constructed_from") <- "formula"
  }else{
    stop(paste0("design argment of class ", class(design), " is not supported. Please ",
                "specify a `model_matrix`, a `character vector`, or a `formula`."))
  }

  if(nrow(model_matrix) != ncol(data)) stop("Number of rows in col_data does not match number of columns of data.\n",
                                            "Were there maybe 'NA's in the colData?")
  if(! is.null(rownames(model_matrix)) &&
     ! all(rownames(model_matrix) == as.character(seq_len(nrow(model_matrix)))) && # That's the default rownames
     ! is.null(colnames(data))){
    if(! all(rownames(model_matrix) == colnames(data))){
      if(setequal(rownames(model_matrix), colnames(data))){
        # Rearrange the rows to match the columns of data
        model_matrix <- model_matrix[colnames(data), ,drop=FALSE]
      }else{
        stop("The rownames of the model_matrix / col_data do not match the column names of data.")
      }
    }

  }

  if(ncol(model_matrix) >= ncol(data)){
    stop("The model_matrix:\n", format_matrix(head(model_matrix)), "\nhas more columns (", ncol(model_matrix),
         ") than the there are samples in the data matrix\n",
         format_matrix(head(data[,seq_len(min(5, ncol(data))),drop=FALSE], n = 3)),
         "\nwhich has ", ncol(data), " columns. ",
         "Too few replicates / too many coefficients to fit model.")
  }

  # Check rank of model_matrix
  qr_mm <- qr(model_matrix)
  if(qr_mm$rank < ncol(model_matrix) && n_samples > 0){
    is_zero_column <- DelayedMatrixStats::colCounts(model_matrix, value = 0) == nrow(model_matrix)
    if(any(is_zero_column)){
      stop("The model matrix:\n", format_matrix(head(model_matrix, n = 3)), "\nseems degenerate ('matrix_rank(model_matrix) < ncol(model_matrix)'). ",
           "Column ", paste0(which(is_zero_column), collapse = ", "), " contains only zeros.")
    }else{
      stop("The model matrix:\n", format_matrix(head(model_matrix, n = 3)), "\nseems degenerate ('matrix_rank(model_matrix) < ncol(model_matrix)'). ",
           "Some columns are perfectly collinear. Did you maybe include the same coefficient twice?")
    }
  }

  rownames(model_matrix) <- colnames(data)
  validate_model_matrix(model_matrix, data)
  model_matrix <- add_attr_if_intercept(model_matrix)
  list(model_matrix = model_matrix, design_formula = design_formula)
}

add_attr_if_intercept <- function(model_matrix){
  intercept_position <- 0
  for(col_idx in seq_len(ncol(model_matrix))){
    has_intercept <- all(model_matrix[,col_idx] == 1)
    if(has_intercept){
      intercept_position <- col_idx
      break
    }
  }
  attr(model_matrix, "intercept_position") <- intercept_position
  model_matrix
}



handle_subsample_parameter <- function(data, subsample){
  if(isFALSE(subsample)){
    n_subsamples <- ncol(data)
  }else if(isTRUE(subsample)){
    n_subsamples <- 1000
  }else{
    n_subsamples <- subsample
  }
  min(n_subsamples, ncol(data))
}


handle_ridge_penalty_parameter <- function(ridge_penalty, model_matrix, verbose){
  if(! is.null(ridge_penalty)){
    ridge_penalty[abs(ridge_penalty) < 1e-10] <- 1e-10

    intercept_position <- attr(model_matrix, "intercept_position")
    if(length(ridge_penalty) == 1){
      ridge_penalty <- rep_len(ridge_penalty, ncol(model_matrix))
      if(! is.null(intercept_position) && intercept_position[1] != 0){
        ridge_penalty[intercept_position] <- 1e-10
      }
    }else if(is.matrix(ridge_penalty)){
      if(ncol(ridge_penalty) != ncol(model_matrix)){
        stop("'ridge_penalty' is a matrix, but its dimensions do not match ",
             "the column of the model_matrix (", ncol(model_matrix), ").")
      }

    }else if(length(ridge_penalty) == ncol(model_matrix)){
      # Got a full length ridge_penalty, check if this conflicts with intercept
      if(! is.null(intercept_position) && any(abs(ridge_penalty[intercept_position]) > 1e-10)){
        warning("A ridge penalty for each column of the design matrix was provided, including the intercept ",
                "in column ", intercept_position, ". Are you sure this is correct?\n",
                "To avoid this message, set the ridge_penalty[", intercept_position, "] to a value ",
                "smaller than 1e-10.")
      }
      ridge_penalty
    }else{
      stop("The definition of the ridge penalty does not match the model_matrix. ",
           "It must either be of length 1 or the number of columns in the design matrix.\n",
           "If length(ridge_penalty) == 1, it is applied to all columns except the intercept.")
    }
  }
  ridge_penalty
}




validate_model_matrix <- function(matrix, data){
  stopifnot(is.matrix(matrix))
  stopifnot(nrow(matrix) == ncol(data))
}


relevel_columns_to_reference_level <- function(col_data, reference_level = NULL){
  if(! is.null(reference_level)){
    has_ref_level <- vapply(col_data, function(x){
      any(!is.na(x) & x == reference_level)
    }, FUN.VALUE = FALSE)
    if(all(has_ref_level == FALSE)){
      stop("None of the columns contains the specified reference_level.")
    }
    col_data[has_ref_level] <- lapply(col_data[has_ref_level], function(col){
      if(is.character(col)){
        col <- as.factor(col)
      }
      stats::relevel(col, ref = reference_level)
    })
  }
  col_data
}


convert_chr_vec_to_model_matrix <- function(design, reference_level){
  if(! is.factor(design)){
    design_fct <- as.factor(design)
  }else{
    design_fct <- design
  }

  if(length(levels(design_fct)) == 1){
    helper_df <- data.frame(x_ = design_fct)
    tmp <- convert_formula_to_model_matrix(~ 1, col_data = helper_df)
    mm <- tmp$model_matrix
    formula <- tmp$formula
    colnames(mm) <- levels(design_fct)
  }else if(is.null(reference_level)){
    helper_df <- data.frame(x_ = design_fct)
    tmp <- convert_formula_to_model_matrix(~ x_ - 1, col_data = helper_df)
    mm <- tmp$model_matrix
    formula <- tmp$formula
    colnames(mm) <- sub("^x_", "", colnames(mm))
  }else{
    design_fct <- stats::relevel(design_fct, ref = reference_level)
    helper_df <- data.frame(x_ = design_fct)
    tmp <- convert_formula_to_model_matrix(~ x_ + 1, col_data = helper_df)
    mm <- tmp$model_matrix
    formula <- tmp$formula
    colnames(mm)[-1] <- paste0(sub("^x_", "", colnames(mm)[-1]), "_vs_", reference_level)
  }
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  list(formula = formula, model_matrix = mm)
}

convert_formula_to_model_matrix <- function(formula, col_data){
  tryCatch({
    mf <- model.frame(formula, data = col_data, drop.unused.levels = TRUE)
    terms <- attr(mf, "terms")
    attr(terms, "xlevels") <- stats::.getXlevels(terms, mf)
    mm <- stats::model.matrix.default(terms, mf)
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Error while parsing the formula (", formula, ").\n",
           "Variable '", match[2], "' not found in col_data or global environment. Possible variables are:\n",
           paste0(colnames(col_data), collapse = ", "), call. = FALSE)
    }else{
      stop(e$message)
    }
  })

  # Otherwise every copy of the model stores the whole global environment!
  attr(terms, ".Environment") <- c()
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  list(formula = terms, model_matrix = mm)
}

extract_data_from_formula <- function(formula, col_data, encl = parent.frame()){
  if(length(formula) < 3){
    stop("The formula does not have a left hand side. Please call this function like this:\n",
         "'glm_gp(data = mat, design = ~ a + b + c, ...)'", call. = FALSE)
  }
  data <- eval(formula[[2]], envir = col_data, enclos = encl)
  formula[[2]] <- NULL
  list(data = data, design = formula)
}


is_on_disk.glmGamPoi <- function(fit){
  is(fit$Mu, "DelayedMatrix") && is(DelayedArray::seed(fit$Mu), "HDF5ArraySeed")
}


