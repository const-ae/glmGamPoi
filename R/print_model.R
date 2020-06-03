
#' Convert glmGamPoi object to a list
#'
#' @param x an object with class glmGamPoi
#' @param ... not used
#'
#' @return
#' The method returns a list with the following elements:
#' \describe{
#'   \item{`Beta`}{a matrix with dimensions `nrow(data) x n_coefficients` where `n_coefficients` is
#'   based on the `design` argument. It contains the estimated coefficients for each gene.}
#'   \item{`overdispersions`}{a vector with length `nrow(data)`. The overdispersion parameter for each
#'   gene. It describes how much more the counts vary than one would expect according to the Poisson
#'   model.}
#'   \item{`Mu`}{a matrix with the same dimensions as `dim(data)`. If the calculation happened on
#'   disk, than `Mu` is a `HDF5Matrix`. It contains the estimated mean value for each gene and
#'   sample.}
#'   \item{`size_factors`}{a vector with length `ncol(data)`. The size factors are the inferred
#'   correction factors for different sizes of each sample. They are also sometimes called the
#'   exposure factor.}
#'   \item{`model_matrix`}{a matrix with dimensions `ncol(data) x n_coefficients`. It is build based
#'   on the `design` argument.}
#' }
#'
#'
#' @export
as.list.glmGamPoi <- function(x, ...){
  class(x) <- "list"
  x
}


#' Pretty print the result from glm_gp()
#'
#' @param x the glmGamPoi object
#' @param object the glmGamPoi object that is summarized
#' @param ... additional parameters, currently ignored
#'
#' @return
#'   The `print()` methods return the object `x`. The
#'   `format()` method returns a string. The `summary()`
#'   method returns an object of class `summary.glmGamPoi`.
#'
#' @export
print.glmGamPoi <- tools::.print.via.format

#' @rdname print.glmGamPoi
#' @export
format.glmGamPoi <- function(x, ...){
  string_builder <- ""
  nrow <- nrow(x$Mu)
  ncol <- ncol(x$Mu)
  npara <- ncol(x$Beta)
  string_builder <- paste0(string_builder, "glmGamPoiFit object:\n",
                           "The data had ", nrow, " rows and ", ncol, " columns.\n",
                           "A model with ", npara, " coefficient was fitted.")
  string_builder
}

#' @rdname print.glmGamPoi
#' @export
summary.glmGamPoi <- function(object, ...){
  ans <- object
  class(ans) <- "summary.glmGamPoi"
  ans
}

#' @rdname print.glmGamPoi
#' @export
print.summary.glmGamPoi <- tools::.print.via.format

#' @rdname print.glmGamPoi
#' @export
format.summary.glmGamPoi <- function(x, ...){
  header <- paste0(format.glmGamPoi(x), "\n")
  nrow <- nrow(x$Mu)
  ncol <- ncol(x$Mu)
  npara <- ncol(x$Beta)
  ndigits <- 3

  # Make string for design formula
  if(is.null(x$design_formula)){
    design_summary <- ""
  }else{
    design_summary <- paste0("The design formula is: Y", format(x$design_formula), "\n")
  }

  # Make String for Beta
  if(nrow < 5 && npara <= 3){
    beta_est_summary <- paste0("Beta:\n", format_matrix(t(x$Beta), digits = ndigits), "\n")
  }else if(nrow < 5){
    beta_est_summary <- paste0("Beta:\n", format_matrix(t(x$Beta[,seq_len(3)]), digits = ndigits), "\n...\n")
  }else if(npara <= 3){
    beta_est_summary <- paste0("Beta:\n", format_matrix(my_matrix_summary(x$Beta), digits = ndigits), "\n")
  }else{
    beta_est_summary <- paste0("Beta:\n", format_matrix(my_matrix_summary(x$Beta[,seq_len(3)]), digits = ndigits), "\n...\n")
  }

  # Make String for Deviance
  deviance_mat <- matrix(x$deviances, ncol = 1, dimnames = list(NULL, ""))
  dev_summary <- paste0("\ndeviance:\n", format_matrix(my_matrix_summary(deviance_mat), digits = ndigits), "\n")

  # Make String for Overdispersion
  overdisp_mat <- matrix(x$overdispersions, ncol = 1, dimnames = list(NULL, ""))
  overdisp_summary <- paste0("\noverdispersion:\n", format_matrix(my_matrix_summary(overdisp_mat), digits = ndigits), "\n")

  # Make String for Shrunken overdispersion
  if(is.null(x$overdispersion_shrinkage_list)){
    disp_shrunk_summary <- "\nNo quasi-likelihood overdispersion available\n"
  }else{
    disp_shrunk_mat <- matrix(x$overdispersion_shrinkage_list$ql_disp_shrunken, ncol = 1, dimnames = list(NULL, ""))
    disp_shrunk_summary <- paste0("\nShrunken quasi-likelihood overdispersion:\n", format_matrix(my_matrix_summary(disp_shrunk_mat), digits = ndigits), "\n")
  }


  # Make String for size factor
  sf_mat <- matrix(x$size_factors, ncol = 1, dimnames = list(NULL, ""))
  sf_summary <- paste0("\nsize_factors:\n", format_matrix(my_matrix_summary(sf_mat), digits = ndigits), "\n")

  # Make String for size factor
  sf_mat <- matrix(x$size_factors, ncol = 1, dimnames = list(NULL, ""))
  sf_summary <- paste0("\nsize_factors:\n", format_matrix(my_matrix_summary(sf_mat), digits = ndigits), "\n")

  # Make String for Mu
  if(is.matrix(x$Mu)){
    mu_mat <- matrix(c(x$Mu), ncol = 1, dimnames = list(NULL, ""))
    mu_summary <- paste0("\nMu:\n", format_matrix(my_matrix_summary(mu_mat), digits = ndigits), "\n")
  }else{
    mu_summary <- paste0("\nMu is a ", class(x$Mu)[1], ", will skip summary to avoid ",
                         "lengthy computation.")
  }

  paste0(header, design_summary, "\n", beta_est_summary, dev_summary, overdisp_summary,
         disp_shrunk_summary, sf_summary, mu_summary)
}


#' Helper to format a matrix nicely
#'
#' @return a string
#'
#' @keywords internal
format_matrix <- function(matrix, digits = NULL){

  rownames <- if(is.null(rownames(matrix))){
    if(nrow(matrix) > 0){
      paste0("[",seq_len(nrow(matrix)), ",]")
    }else{
      character(0)
    }
  }else{
    rownames(matrix)
  }
  colnames <- if(is.null(colnames(matrix))){
    if(ncol(matrix) > 0){
      paste0("[,",seq_len(ncol(matrix)), "]")
    }else{
      character(0)
    }
  }else{
    colnames(matrix)
  }
  offset <- max(nchar(rownames), 0)
  char_matrix <- vapply(seq_len(length(colnames)), function(col_idx){
    format(matrix[,col_idx], digits = digits)
  }, FUN.VALUE = rep("", nrow(matrix)))
  if(is.vector(char_matrix)){
    char_matrix <- matrix(char_matrix, nrow = 1)
  }
  width_per_col <- vapply(seq_len(length(colnames)), function(col_idx){
    max(nchar(char_matrix[, col_idx]), 0)
  }, FUN.VALUE = 0.0)
  width_per_header <- nchar(colnames)
  width_per_col <- pmax(width_per_col, width_per_header)
  adapted_colnames <- vapply(seq_len(length(colnames)), function(col_idx){
    format(colnames[col_idx], width = width_per_col[col_idx], justify = "right")
  }, FUN.VALUE = "")
  head_line <- paste0(c(format("", width = offset), adapted_colnames), collapse = " ")
  rows_pasted <- vapply(seq_len(length(rownames)), function(row_idx){
    paste0(vapply(seq_len(length(colnames)), function(col_idx){
      format(char_matrix[row_idx, col_idx], width = width_per_col[col_idx], justify = "right")
    }, FUN.VALUE = ""), collapse = " ")
  }, FUN.VALUE = "")
  adapted_rownames <- vapply(seq_len(length(rownames)), function(row_idx){
    format(rownames[row_idx], width = offset, justify = "right")
  }, FUN.VALUE = "")
  paste0(head_line, "\n", paste0(adapted_rownames, " ", rows_pasted, collapse = "\n"))
}



my_vec_summary <- function(x, count_na = TRUE){
  res <- unname(quantile(x, probs = c(0, 0.25, 0.5, .75, 1), na.rm=TRUE))
  names(res) <- c("Min", "1st Qu.", "Median", "3rd Qu.", "Max")
  if(count_na){
    res <- c(res, `#NA` = sum(is.na(x)))
  }
  res
}

my_matrix_summary <- function(X){
  if(any(is.na(X))){
    res <- t(vapply(seq_len(ncol(X)), function(col_idx){
      my_vec_summary(X[, col_idx], count_na = TRUE)
    }, FUN.VALUE = rep(0.0, 6)))
  }else{
    res <- t(vapply(seq_len(ncol(X)), function(col_idx){
      my_vec_summary(X[, col_idx], count_na = FALSE)
    }, FUN.VALUE = rep(0.0, 5)))
  }
  rownames(res) <- colnames(X)
  res
}
