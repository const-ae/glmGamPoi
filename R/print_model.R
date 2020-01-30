
print.glmGamPoi <- tools::.print.via.format

#' Pretty print the result from glm_gp()
#'
#' @keywords internal
format.glmGamPoi <- function(x, ...){
  string_builder <- ""
  nrow <- nrow(x$Mu_est)
  ncol <- ncol(x$Mu_est)
  npara <- ncol(x$Beta_est)
  string_builder <- paste0(string_builder, "glmGamPoiFit object:\n")
  string_builder <- paste0(string_builder, "The data had ", nrow, " rows and ", ncol,
                           " columns.\n")
  string_builder <- paste0(string_builder, "A model with ", npara, " coefficient was fitted.")
  string_builder
}

summary.glmGamPoi <- function(object, ...){
  ans <- object
  class(ans) <- "summary.glmGamPoi"
  ans
}

print.summary.glmGamPoi <- tools::.print.via.format

format.summary.glmGamPoi <- function(x, ...){
  header <- paste0(format.glmGamPoi(x), "\n")
  nrow <- nrow(x$Mu_est)
  ncol <- ncol(x$Mu_est)
  npara <- ncol(x$Beta_est)
  ndigits <- 3

  # Make string for design formula
  if(is.null(x$design_formula)){
    design_summary <- ""
  }else{
    design_summary <- paste0("The design formula is: Y", format(x$design_formula), "\n")
  }

  # Make String for Beta_est
  if(nrow < 5 && npara <= 3){
    beta_est_summary <- paste0("Beta_est:\n", format_matrix(t(x$Beta_est), digits = ndigits), "\n")
  }else if(nrow < 5){
    beta_est_summary <- paste0("Beta_est:\n", format_matrix(t(x$Beta_est[,1:3]), digits = ndigits), "\n...\n")
  }else if(npara <= 3){
    beta_est_summary <- paste0("Beta_est:\n", format_matrix(my_matrix_summary(x$Beta_est), digits = ndigits), "\n")
  }else{
    beta_est_summary <- paste0("Beta_est:\n", format_matrix(my_matrix_summary(x$Beta_est[,1:3]), digits = ndigits), "\n...\n")
  }

  # Make String for Overdispersion
  overdisp_mat <- matrix(x$overdispersions, ncol = 1, dimnames = list(NULL, ""))
  overdisp_summary <- paste0("\noverdispersion:\n", format_matrix(my_matrix_summary(overdisp_mat), digits = ndigits), "\n")

  # Make String for size factor
  sf_mat <- matrix(x$size_factors, ncol = 1, dimnames = list(NULL, ""))
  sf_summary <- paste0("\nsize_factors:\n", format_matrix(my_matrix_summary(sf_mat), digits = ndigits), "\n")

  # Make String for size factor
  sf_mat <- matrix(x$size_factors, ncol = 1, dimnames = list(NULL, ""))
  sf_summary <- paste0("\nsize_factors:\n", format_matrix(my_matrix_summary(sf_mat), digits = ndigits), "\n")

  # Make String for Mu_est
  mu_mat <- matrix(c(x$Mu_est), ncol = 1, dimnames = list(NULL, ""))
  mu_summary <- paste0("\nMu_est:\n", format_matrix(my_matrix_summary(mu_mat), digits = ndigits), "\n")


  paste0(header, design_summary, "\n", beta_est_summary, overdisp_summary, sf_summary, mu_summary)
}

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