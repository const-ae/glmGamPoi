



#' Estimate the Size Factors
#'
#' @param Y any matrix-like object (\code{base::matrix()}, \code{DelayedArray}, \code{HDF5Matrix},
#'   \code{Matrix::Matrix()}, etc.)
#'
#'
#' @return a vector with one size factor per column of `Y`
#'
#' @keywords internal
estimate_size_factors <- function(Y){
  # Accept any matrix-like object
  log_geometric_means <- DelayedMatrixStats::rowMeans2(log(Y + 0.5))
  Y2 <- Y
  Y2[Y2 == 0] <- NA
  sf <- exp(DelayedMatrixStats::colMedians(subtract_vector_from_each_column(log(Y2), log_geometric_means), na.rm = TRUE))
  all_zero_column <- is.nan(sf)
  if(any(all_zero_column)){
    sf <- sf/exp(mean(log(sf), na.rm=TRUE))
    sf[all_zero_column] <- 0.001
  }else{
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}



combine_size_factors_and_offset <- function(offset, size_factors, Y, verbose = FALSE){
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)

  make_offset_hdf5_mat <- is(Y, "DelayedMatrix") && is(DelayedArray::seed(Y), "HDF5ArraySeed")

  if(is.matrix(offset)){
    stopifnot(dim(offset) == c(n_genes, n_samples))
    offset_matrix <- offset
  }else{
    stopifnot(length(offset) == 1 || length(offset) == n_samples)
    if(make_offset_hdf5_mat){
      offset_matrix <- DelayedArray::DelayedArray(DelayedArray::SparseArraySeed(c(n_genes, n_samples)))
      offset_matrix <- add_vector_to_each_row(offset_matrix, offset)
    }else{
      offset_matrix <- matrix(offset, nrow=n_genes, ncol = n_samples, byrow = TRUE)
    }
  }
  if(isTRUE(size_factors)){
    if(n_genes <= 1){
      warning("Calculating the size factor does not make sense for a single gene.\n",
              "Consider setting `size_factors = FALSE`.\n",
              "Calculation will however proceed, be careful when interpreting the coefficients.",
              call. = FALSE)
    }
    if(verbose){ message("Calculate Size Factors") }
    lsf <- log(estimate_size_factors(Y))
  }else if(isFALSE(size_factors)){
    lsf <- rep(0, n_samples)
  }else{
    stopifnot(is.numeric(size_factors) && (length(size_factors) == 1 || length(size_factors) == n_samples))
    if(any(size_factors < 0)){
      stop("size factor 'size_factors' must be larger than 0")
    }
    if(length(size_factors) == 1){
      lsf <- rep(log(size_factors), n_samples)
    }else{
      lsf <- log(size_factors)
    }
  }
  # offset_matrix <- offset_matrix + lsf_mat
  offset_matrix <- add_vector_to_each_row(offset_matrix, lsf)
  if(make_offset_hdf5_mat){
    offset_matrix <- HDF5Array::writeHDF5Array(offset_matrix)
  }
  list(offset_matrix = offset_matrix, size_factors = exp(lsf))
}


