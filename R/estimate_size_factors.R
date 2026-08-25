



#' Estimate the Size Factors
#'
#' @param Y any matrix-like object (\code{base::matrix()}, \code{DelayedArray}, \code{HDF5Matrix},
#'   \code{Matrix::Matrix()}, etc.)
#' @param method one of \code{c("normed_sum", "deconvolution", "poscounts")}
#'
#'
#' @return a vector with one size factor per column of `Y`
#'
#' @keywords internal
estimate_size_factors <- function(Y, method, verbose = FALSE){
  if(nrow(Y) <= 1){
    if(verbose) {
      message("Skipping size factor estimation, because there is only a single gene.",
              call. = FALSE)
    }
    return(rep(1, ncol(Y)))
  }
  stopifnot(length(method) == 1 && is.character(method))

  if(method == "poscounts"){
    # Accept any matrix-like object
    log_geometric_means <- DelayedMatrixStats::rowMeans2(log(Y + 0.5))
    Y2 <- Y
    Y2[Y2 == 0] <- NA
    sf <- exp(DelayedMatrixStats::colMedians(subtract_vector_from_each_column(log(Y2), log_geometric_means), na.rm = TRUE))
  }else if(method == "deconvolution"){
    if(requireNamespace("scran", quietly = TRUE)){
      tryCatch({
        sf <- scran::calculateSumFactors(Y)
      }, error = function(err){
        stop("Error in size factor estimation with 'size_factors = \"deconvolution\"'.\n",
             "'scran::calculateSumFactors(Y)' threw the following error: \n\t",
             err, call. = FALSE)
      })
    }else{
      stop("To use the \"deconvolution\" method for size factor calculation, you need to install the ",
           "'scran' package from Bioconductor. Alternatively, you can use \"normed_sum\" or \"poscounts\"",
           "to calculate the size factors.", call. = FALSE)
    }
  }else if(method == "normed_sum"){
    sf <- DelayedMatrixStats::colSums2(Y)
  }else if(method == "ratio"){
    log_geo_means <- DelayedMatrixStats::rowMeans2(log(Y))
    sf <- apply(Y, 2, function(cnts) {
      exp(median((log(cnts) - log_geo_means)[is.finite(log_geo_means) & cnts > 0]))
    })
  }else{
    stop("Unknown size factor estimation method: ", method)
  }



  # stabilize size factors to have geometric mean of 1
  all_zero_column <- is.na(sf) | sf <= 0
  sf[all_zero_column] <- NA
  if(any(all_zero_column) && method == "ratio"){
    stop("Too many zeros to calculate the size factors. Please use 'size_factors = \"normed_sum\" or 'size_factors = \"deconvolution\".")
  }
  if(any(all_zero_column)){
    warning(sum(all_zero_column), " columns contain too many zeros to calculate a size factor. The size factor will be fixed to 0.001")
    sf <- sf/exp(mean(log(sf), na.rm=TRUE))
    sf[all_zero_column] <- 0.001
  }else{
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}



combine_size_factors_and_offset <- function(offset, size_factors, Y, verbose = FALSE, offset_as_vec = FALSE){
  n_genes <- nrow(Y)
  n_samples <- ncol(Y)

  zero_offset <- TRUE
  make_offset_hdf5_mat <- is(Y, "DelayedMatrix") && is(DelayedArray::seed(Y), "HDF5ArraySeed")
  make_offset_delayed_mat <- is(Y, "DelayedMatrix")

  if(is.numeric(offset) && ! is.matrix(offset) && (make_offset_hdf5_mat || is(Y, "DelayedMatrix"))){
    zero_offset <- all(offset == 0)
    offset_matrix <-if(length(offset) == 1){
      DelayedArray::ConstantArray(c(n_genes, n_samples), value = offset)
    }else if(length(offset) == n_samples){
      DelayedArray::RleArray(S4Vectors::Rle(offset, lengths = n_genes), dim = c(n_genes, n_samples))
    }else{
      stop("'offset' is a vector. Its length must either be 1 or match the number of samples (", n_samples, ").")
    }
  }else if(is.numeric(offset) && ! is.matrix(offset)){
    if(length(offset) != 1 && length(offset) != n_samples){
      stop("'offset' is a vector. Its length must either be 1 or match the number of samples (", n_samples, ").")
    }
    if(offset_as_vec){
      offset_matrix <- NULL
    }else{
      offset_matrix <- matrix(offset, nrow = n_genes, ncol = n_samples, byrow = TRUE)
    }
  }else if(is.matrix(offset) || is(offset, "DelayedMatrix")){
    stopifnot(dim(offset) == c(n_genes, n_samples))
    offset_matrix <- offset
  }else{
    stop("'offset' is of type ", class(offset)[1], ", which is not supported.")
  }


  if(isTRUE(size_factors) || is.character(size_factors)){
    if(! zero_offset){
      warning("The offset is non zero, nonetheless an additional size factor is estimated. The effective 'offset' is thus 'offset + log(size_factor)'.\n",
              "To suppress this warning either set 'size_factor = FALSE' or 'offset = 0'.")
    }
    method <- if(isTRUE(size_factors)){
      if(requireNamespace("scran", quietly = TRUE)){
        "deconvolution"
      }else{
        "normed_sum"
      }
    }else if(all(size_factors == c("normed_sum", "deconvolution", "poscounts", "ratio"))){
      "normed_sum"
    }else if(length(size_factors) == 1 && size_factors %in% c("normed_sum", "deconvolution", "poscounts", "ratio")){
      size_factors
    }else{
      stop("Cannot handle size factor ", paste0(size_factors, collapse = ", "))
    }
    if(verbose){ message("Calculate Size Factors (", method, ")") }
    lsf <- log(estimate_size_factors(Y, method = method, verbose = verbose))
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
  if(is.null(offset_matrix)){
    offset_val <- offset + lsf
  }else{
    offset_val <- add_vector_to_each_row(offset_matrix, lsf)
  }

  if(make_offset_hdf5_mat && ! is(offset_val, "HDF5Array")){
    # @TODO: figure out a better way to construct a HDF5Array in this case?
    if(is.null(offset_matrix)){
      offset_val <- matrix(offset_val, nrow = n_genes, ncol = n_samples, byrow=TRUE)
    }
    offset_val <- HDF5Array::writeHDF5Array(offset_val)
  }else if(make_offset_delayed_mat && ! is(offset_val, "DelayedArray")){
    # @TODO: figure out a better way to construct a DelayedArray in this case?
    if(is.null(offset_matrix)){
      offset_val <- matrix(offset_val, nrow = n_genes, ncol = n_samples, byrow=TRUE)
    }
    offset_val <- DelayedArray::DelayedArray(offset_val)
  }

  list(offset_matrix = offset_val, size_factors = exp(lsf))
}


