
#' Create a 'SingleCellExperiment' containing pseudo-bulk samples
#'
#' @param sce a 'SingleCellExperiment' or an object of a related class
#' @param group_by an unquoted expression that can refer to columns in
#'   the 'colData()'. All observations with the same factor level are aggregated.
#'   The argument follows the same logic as `dplyr::group_by()`. The argument must
#'   wrapped using `vars()`.
#' @param ... named expressions that summarize columns in 'colData()'. Each expression
#'   must produce a value of length 1. The arguments follow the same logic
#'   as `dplyr::summarize()`.
#' @param aggregation_functions a named list with functions that are used to
#'   aggregate the assays in the `sce`.
#'
#' @return a SingleCellExperiment object
#'
#' @examples
#'  library(SingleCellExperiment)
#'  data <- data.frame(sample = sample(c("sample1", "sample2", "sample3"), size = 50, replace = TRUE),
#'                     celltype = sample(c("T cells", "B cells", "Macrophages"), size = 50, replace = TRUE),
#'                     size = rnorm(n = 50, mean = 40, sd = 15))
#'  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
#'  sce <- SingleCellExperiment(Y, colData = data)
#'  aggr_sce <- pseudobulk_sce(sce, group_by = vars(sample, celltype), size = mean(size))
#'  aggr_sce
#'  colData(aggr_sce)
#'
#' @export
pseudobulk_sce <- function(sce, group_by, ...,
                           aggregation_functions = list(counts = "rowSums2", .default = "rowMeans2")){
  if(missing(group_by)) stop("The 'group_by' must not be missing.")
  tryCatch({
    if(inherits(group_by, "uneval") || !rlang::is_quosures(group_by)) stop("The 'group_by' argument must be wrapped using 'vars()'")
  }, error = function(e){
    stop(e, "\n", "The 'group_by' argument must be wrapped using 'vars()'")
  })

  index_seq <- seq_len(ncol(sce))
  col_data <- as.data.frame(SummarizedExperiment::colData(sce))
  if(! ".default" %in% names(aggregation_functions)){
    aggregation_functions$.default <- "rowMeans"
  }

  # Make vector that is used to
  groups <- lapply(group_by, rlang::eval_tidy, data = col_data)
  if(! all(lengths(groups) == 1 | lengths(groups) == ncol(sce))){
    stop("The argument 'group_by' has lengths ", paste0(lengths(groups), collapse = ","), ", which does not match the number of columns ",
         "in 'sce' (", ncol(sce), ")")
  }else if(any(lengths(groups) == 1)){
    groups[lengths(groups) == 1] <- lapply(groups[lengths(groups) == 1], \(x) rep_len(x, ncol(sce)))
  }

  if(is.null(groups)){
    stop("'group_by' must not be 'NULL'.")
  }else{
    group_split <- split(index_seq, groups, drop = TRUE)
  }

  # Aggregate all assays
  new_assays <- lapply(SummarizedExperiment::assayNames(sce), function(assay_name){
    aggr_fnc <- get_aggregation_function(assay_name, aggregation_functions)
    data_mat <- SummarizedExperiment::assay(sce, assay_name)
    new_data_mat <- do.call(cbind, lapply(group_split, function(idx){
      aggr_fnc(data_mat[,idx,drop=FALSE])
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
      new_data_mat <- do.call(cbind, lapply(group_split, function(idx){
        aggr_fnc(data_mat[,idx,drop=FALSE])
      }))
      SingleCellExperiment::LinearEmbeddingMatrix(t(new_data_mat), SingleCellExperiment::featureLoadings(tdata_mat),
                                                  factorData = SingleCellExperiment::factorData(tdata_mat))
    }else{
      data_mat <- t(tdata_mat)
      new_data_mat <- do.call(cbind, lapply(group_split, function(idx){
        aggr_fnc(data_mat[,idx,drop=FALSE])
      }))
      rownames(new_data_mat) <- rownames(data_mat)
      t(new_data_mat)
    }
  })
  names(new_red_dims) <- SingleCellExperiment::reducedDimNames(sce)

  # Aggregate column data
  dots_cap <- rlang::enquos(...)
  dots <- lapply(dots_cap, rlang::eval_tidy, data = col_data)
  new_col_data <- lapply(seq_along(dots), function(dot_idx){
    dot <- dots[[dot_idx]]
    dot_name <- names(dots_cap)[dot_idx]
    if(is.null(dot_name) || dot_name == ""){
      rlang::as_label(dot)
    }
    unname(do.call(c, lapply(group_split, function(idx){
      res <- rlang::eval_tidy(dot, col_data[idx,,drop=FALSE])
      if(length(res) != 1) stop("Illegal result in aggregation of '", dot_name,
                                "'. The aggregated value has to be of length 1. However, it has length ", length(res))
      res
    })))
  })
  if(is.null(names(dots_cap))){
    names(new_col_data) <- vapply(dots_cap, rlang::as_label, FUN.VALUE = character(1L))
  }else{
    names(new_col_data) <- names(dots_cap)
  }
  names(new_col_data)[names(new_col_data) == ""] <- vapply(dots_cap[names(new_col_data) == ""], rlang::as_label, FUN.VALUE = character(1L))

  # Make id columns
  row_sel <- vapply(group_split, head, n = 1, FUN.VALUE = integer(1L))
  id_columns <- as.data.frame(groups)[row_sel,,drop=FALSE]
  id_columns <- lapply(seq_along(groups), \(idx){
    factor(id_columns[[idx]], levels = levels(as.factor(groups[[idx]])))
  })
  names(id_columns) <-  vapply(seq_along(group_by), \(idx){
    id_name <- names(group_by)[idx]
    if(is.null(id_name) || id_name == ""){
      id_name <- rlang::as_label(group_by[[idx]])
    }
    id_name
  }, FUN.VALUE = character(1L))
  new_col_data <- c(id_columns, new_col_data)

  SingleCellExperiment::SingleCellExperiment(new_assays, colData = new_col_data, reducedDims = new_red_dims)
}

get_aggregation_function <- function(assay_name, aggregation_functions){
  aggr_fnc <- if(assay_name %in% names(aggregation_functions)){
    aggregation_functions[[assay_name]]
  }else{
    aggregation_functions[[".default"]]
  }
  if(is.character(aggr_fnc)){
    aggr_fnc <- if(aggr_fnc == "rowSums2"){
      MatrixGenerics::rowSums2
    }else if(aggr_fnc == "rowMeans2"){
      MatrixGenerics::rowMeans2
    }else{
      get(aggr_fnc, envir =  globalenv(), mode = "function")
    }
  }

  aggr_fnc
}

#' Quote grouping variables
#'
#' @param ... the quoted expression
#'
#' @seealso ggplot2::vars, dplyr::vars
#'
#' @export
vars <- function(...){
  rlang::quos(...)
}

