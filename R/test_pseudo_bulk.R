

test_pseudo_bulk <- function(data, design,
                             aggregate_cells_by,
                             contrast,
                             reduced = NULL,
                             subset_to = NULL, col_data = NULL, reference_level = NULL,
                             pval_adjust_method = "BH", sort_by = NULL,
                             decreasing = FALSE, n_max = Inf,
                             verbose = FALSE){
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
  pseudo_bulk_split <- split(seq_len(ncol(data))[subset_to], col_data[[aggregate_cells_by]][subset_to])

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



  fit <- glm_gp(new_data_mat,
                design = new_model_matrix,
                verbose = verbose)

  cnt_capture <- substitute(contrast)
  test_de(fit, contrast = cnt_capture, reduced = reduced,
          subset = NULL, pseudo_bulk = NULL,
          pval_adjust_method =pval_adjust_method, sort_by = sort_by,
          decreasing = decreasing, n_max = n_max,
          verbose = verbose)

}
