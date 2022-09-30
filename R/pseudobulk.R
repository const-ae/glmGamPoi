
#' Create a 'SingleCellExperiment' containing pseudo-bulk samples
#'
#' @param sce a 'SingleCellExperiment' or an object of a related class
#' @param pseudobulk_by an unquoted expression that can refer to columns in
#'   the 'colData()'. All observations with the same factor level are aggregated.
#'   The argument follows the same logic as `dplyr::group_by()`.
#' @param ... named expressions that summarize columns in 'colData()'. Each expression
#'   must produce a value of length 1. The arguments follow the same logic
#'   as `dplyr::summarize()`.
#' @param aggregation_functions a named list with functions that are used to
#'   aggregate the assays in the `sce`.
#'
#'
#' @export
pseudobulk_sce <- function(sce, pseudobulk_by, ...,
                           aggregation_functions = list(counts = "rowSums", .default = "rowMeans")){
  if(missing(pseudobulk_by)) stop("The 'pseudobulk_by' must not be missing.")
  pseudobulk_by_capture <- substitute(pseudobulk_by)
  pseudobulk_sce_q(sce, pseudobulk_by = pseudobulk_by_capture, ...,
                   aggregation_functions = aggregation_functions, env = parent.frame())
}

pseudobulk_sce_q <- function(sce, pseudobulk_by, ..., aggregation_functions, env = parent.frame()){

  index_seq <- seq_len(ncol(sce))
  col_data <- SummarizedExperiment::colData(sce)
  if(! ".default" %in% names(aggregation_functions)){
    aggregation_functions$.default <- "rowMeans"
  }

  # Make vector that is used to
  pseudobulk_by_e <- eval_with_q(pseudobulk_by, col_data, env = env)
  if(length(pseudobulk_by_e) != 1 && length(pseudobulk_by_e) != ncol(sce)){
    stop("'pseudobulk_by' is length ", length(pseudobulk_by_e), " but it has to be the same ",
         "length as the number of columns in 'sce' (", ncol(sce), ")")
  }else if(length(pseudobulk_by_e) == 1){
    pseudobulk_by_e <- rep_len(pseudobulk_by_e, ncol(sce))
  }

  if(is.null(pseudobulk_by_e)){
    stop("'pseudobulk_by' must not be 'NULL'.")
  }else{
    pseudo_bulk_split <- split(index_seq, pseudobulk_by_e, drop = TRUE)
  }

  # Aggregate all assays
  new_assays <- lapply(SummarizedExperiment::assayNames(sce), function(assay_name){
    aggr_fnc <- get_aggregation_function(assay_name, aggregation_functions)
    data_mat <- SummarizedExperiment::assay(sce, assay_name)
    new_data_mat <- do.call(cbind, lapply(pseudo_bulk_split, function(idx){
      aggr_fnc(data_mat, cols = idx)
    }))
    rownames(new_data_mat) <- rownames(sce)
    new_data_mat
  })
  names(new_assays) <- SummarizedExperiment::assayNames(sce)

  # Aggregate reduced dims
  new_red_dims <- lapply(SingleCellExperiment::reducedDimNames(sce), function(red_name){
    aggr_fnc <- get_aggregation_function(red_name, aggregation_functions)
    tdata_mat <- SingleCellExperiment::reducedDim(sce, red_name)
    if(is(tdata_mat, "LinearEmbeddingMatrix")){
      data_mat <- t(SingleCellExperiment::sampleFactors(tdata_mat))
      new_data_mat <- do.call(cbind, lapply(pseudo_bulk_split, function(idx){
        aggr_fnc(data_mat, cols = idx)
      }))
      SingleCellExperiment::LinearEmbeddingMatrix(t(new_data_mat), SingleCellExperiment::featureLoadings(tdata_mat),
                                                  factorData = SingleCellExperiment::factorData(tdata_mat))
    }else{
      data_mat <- t(tdata_mat)
      new_data_mat <- do.call(cbind, lapply(pseudo_bulk_split, function(idx){
        aggr_fnc(data_mat, cols = idx)
      }))
      rownames(new_data_mat) <- rownames(data_mat)
      t(new_data_mat)
    }
  })
  names(new_red_dims) <- SingleCellExperiment::reducedDimNames(sce)

  # Aggregate column data
  dots <- eval(substitute(alist(...)))
  new_col_data <- lapply(seq_along(dots), function(dot_idx){
    dot <- dots[[dot_idx]]
    dot_name <- if(is.null(names(dots))){
      dot_idx
    }else{
      names(dots)[dot_idx]
    }
    unname(do.call(c, lapply(pseudo_bulk_split, function(idx){
      res <- eval_with_q(dot, col_data[idx,], env = env)
      if(length(res) != 1) stop("Illegal result in aggregation of '", dot_name,
                                "'. The aggregated value has to be of length 1.")
      res
    })))
  })

  if(is.null(names(dots))){
    names(new_col_data) <- vapply(dots, deparse_one, FUN.VALUE = character(1L))
  }else{
    names(new_col_data) <- names(dots)
  }
  names(new_col_data)[names(new_col_data) == ""] <- vapply(dots[names(new_col_data) == ""], deparse_one, FUN.VALUE = character(1L))
  id_column <- list(factor(levels(droplevels(as.factor(pseudobulk_by_e))),
                           levels = levels(as.factor(pseudobulk_by_e))))
  names(id_column) <- deparse_one(pseudobulk_by)
  new_col_data <- c(id_column, new_col_data)

  SingleCellExperiment::SingleCellExperiment(new_assays, colData = new_col_data, reducedDims = new_red_dims)

}

get_aggregation_function <- function(assay_name, aggregation_functions){
  aggr_fnc <- if(assay_name %in% names(aggregation_functions)){
    aggregation_functions[[assay_name]]
  }else{
    aggregation_functions[[".default"]]
  }
  if(is.character(aggr_fnc)){
    aggr_fnc <- if(aggr_fnc == "rowSums"){
      MatrixGenerics::rowSums2
    }else if(aggr_fnc == "rowMeans"){
      MatrixGenerics::rowMeans2
    }else{
      get(aggr_fnc, envir = global_env(), mode = "function")
    }
  }

  aggr_fnc
}


deparse_one <- function(expr){
  str <- deparse(expr, 60L)
  if (length(str) > 1) {
    if (is_call(expr)) {
      str <- deparse(call2(expr[[1]], quote(...)), 60L)
    }
    str <- paste(str, collapse = "\n")
  }
  str
}
