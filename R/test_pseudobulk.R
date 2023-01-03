

test_pseudobulk <- function(data, design,
                            aggregate_cells_by,
                            contrast,
                            reduced_design = NULL,
                            ridge_penalty = 0,
                            subset_to = NULL, col_data = NULL, reference_level = NULL,
                            pval_adjust_method = "BH", sort_by = NULL,
                            decreasing = FALSE, n_max = Inf,
                            verbose = FALSE){

  aggregate_cells_by_capture <- substitute(aggregate_cells_by)
  subset_to_capture <- substitute(subset_to)
  sort_by_capture <- substitute(sort_by)

  test_pseudobulk_q(data, design = design,
                    aggregate_cells_by = aggregate_cells_by_capture,
                    contrast = {{contrast}},
                    reduced_design = reduced_design,
                    ridge_penalty = ridge_penalty,
                    subset_to = subset_to_capture,
                    col_data = col_data,
                    reference_level = reference_level,
                    pval_adjust_method = pval_adjust_method,
                    sort_by = sort_by_capture,
                    decreasing = decreasing, n_max = n_max,
                    verbose = verbose,
                    env = parent.frame())

}

test_pseudobulk_q <- function(data, design,
                              aggregate_cells_by,
                              contrast,
                              reduced_design = NULL,
                              ridge_penalty = 0,
                              subset_to = NULL, col_data = NULL, reference_level = NULL,
                              pval_adjust_method = "BH", sort_by = NULL,
                              decreasing = FALSE, n_max = Inf,
                              verbose = FALSE,
                              env = parent.frame()){


  # Make sure 'subset_to' is valid
  if(is.null(subset_to)){
    subset_to <- TRUE
  }
  # Validate `data`
  if(is.vector(data)){
    data <- matrix(data, nrow = 1)
  }
  data_mat <- handle_data_parameter(data, on_disk = NULL)

  # Convert the formula to a model_matrix
  col_data <- get_col_data(data, col_data)
  des <- handle_design_parameter(design, data, col_data, reference_level)

  # Split the data according to the `aggregate_cells_by` column
  index_seq <- seq_len(ncol(data))
  names(index_seq) <- colnames(data)
  aggregate_cells_by_e <- eval_with_q(aggregate_cells_by, col_data, env = env)
  if(length(aggregate_cells_by_e) == length(ncol(data)) ){
    stop("'aggregate_cells_by' must be exactly as long as the number of columns in data")
  }
  names(aggregate_cells_by_e) <- colnames(data)
  subset_to_e <- eval_with_q(subset_to, col_data, env = env)
  pseudo_bulk_split <- split(index_seq[subset_to_e], aggregate_cells_by_e[subset_to_e])

  # Aggregate the model matrix
  new_model_matrix <- do.call(rbind, lapply(pseudo_bulk_split, function(idx){
    DelayedMatrixStats::colMeans2(des$model_matrix, rows = idx)
  }))
  colnames(new_model_matrix) <- colnames(des$model_matrix)

  # Aggregate the count matrix
  new_data_mat <- do.call(cbind, lapply(pseudo_bulk_split, function(idx){
    DelayedMatrixStats::rowSums2(data_mat, cols = idx)
  }))
  rownames(new_data_mat) <- rownames(data_mat)

  # Aggregate reduced design if not NULL
  if(! is.null(reduced_design)) {
    reduced_mat <- handle_design_parameter(reduced_design, data, col_data, reference_level)$model_matrix
    new_reduced_mat <- do.call(rbind, lapply(pseudo_bulk_split, function(idx){
      DelayedMatrixStats::colMeans2(reduced_mat, rows = idx)
    }))
    colnames(new_reduced_mat) <- colnames(reduced_mat)
    reduced_design <- new_reduced_mat
  }


  fit <- glm_gp(new_data_mat,
                design = new_model_matrix,
                ridge_penalty = ridge_penalty,
                verbose = verbose)

  test_de_q(fit, contrast = {{contrast}}, reduced_design = reduced_design,
            subset_to = NULL, pseudobulk_by = NULL,
            pval_adjust_method = pval_adjust_method, sort_by = sort_by,
            decreasing = decreasing, n_max = n_max,
            verbose = verbose,
            env = env)

}
