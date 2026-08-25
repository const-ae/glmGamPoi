calculate_mu <- function(Beta, model_matrix, offset_matrix) {
  make_offset_hdf5_mat <- is(offset_matrix, "DelayedMatrix")
  if (make_offset_hdf5_mat) {
    mu <- exp(delayed_matrix_multiply(DelayedArray::DelayedArray(Beta), DelayedArray::DelayedArray(t(model_matrix))) + offset_matrix)
    mu <- HDF5Array::writeHDF5Array(mu)
    mu
  } else {
    if (is.matrix(offset_matrix)) {
      exp(Matrix::tcrossprod(Beta, model_matrix) + offset_matrix)
    } else {
      stopifnot(is.vector(offset_matrix, mode = "numeric"))
      exp(add_vector_to_each_row(Matrix::tcrossprod(Beta, model_matrix), offset_matrix))
    }
  }
}

mk_delayed_mu <- function(Beta, model_matrix, offset) {
  fn <- function(i = NULL, j = NULL) {
    if (!is.null(i)) {
      i <- handle_sub_param(nrow(Beta), i)
      Beta <- Beta[i, , drop = FALSE]
      if (!is.vector(offset)) { offset <- offset[i, , drop = FALSE] }
    }
    if (!is.null(j)) {
      j <- handle_sub_param(nrow(model_matrix), j)
      model_matrix <- model_matrix[j, , drop = FALSE]
      offset <- if (is.vector(offset)) { offset[j] } else { offset[, j, drop = FALSE] }
    }
    calculate_mu(Beta, model_matrix, offset)
  }
  attr(fn, "betas") <- Beta
  attr(fn, "model_matrix") <- model_matrix
  attr(fn, "offsets") <- offset

  fn
}
