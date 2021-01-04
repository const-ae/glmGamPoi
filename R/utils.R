
delayed_matrix_apply_block <- function(Y, Mu, overdispersion, FUN){
  stopifnot(nrow(Y) == nrow(Mu), ncol(Y) == ncol(Mu))
  res_sink <- HDF5Array::HDF5RealizationSink(dim(Y))
  on.exit({
    DelayedArray::close(res_sink)
  }, add = TRUE)

  res_grid <- DelayedArray::defaultAutoGrid(res_sink)

  for (coord1 in seq_len(ncol(res_grid))) {
    for(coord2 in seq_len(nrow(res_grid))){
      sel <- res_grid[[coord2, coord1]]
      Y_block <- DelayedArray::read_block(Y, sel)
      Mu_block <- DelayedArray::read_block(Mu, sel)
      res_block <- FUN(Y_block, Mu_block, overdispersion[seq(sel@ranges@width[1]) - 1 + sel@ranges@start[1]])
      DelayedArray::write_block(res_sink, res_grid[[coord2, coord1]], res_block)
    }
  }

  as(res_sink, "DelayedArray")
}



delayed_matrix_multiply <- function(x, y){
  res_sink <- HDF5Array::HDF5RealizationSink(c(nrow(x), ncol(y)))
  on.exit({
    DelayedArray::close(res_sink)
  }, add = TRUE)

  res_grid <- DelayedArray::defaultAutoGrid(res_sink)

  row_ticks <- cumsum(vapply(seq_len(dim(res_grid)[1]), function(idx){
    dim(res_grid[[idx, 1L]])[1]
  }, FUN.VALUE = 0L))
  col_ticks <- cumsum(vapply(seq_len(dim(res_grid)[2]), function(idx){
    dim(res_grid[[1L, idx]])[2]
  }, FUN.VALUE = 0L))



  x_grid <- DelayedArray::ArbitraryArrayGrid(tickmarks = list(row_ticks, ncol(x)))
  y_grid <- DelayedArray::ArbitraryArrayGrid(tickmarks = list(nrow(y), col_ticks))


  for (coord1 in seq_len(ncol(res_grid))) {
    for(coord2 in seq_len(nrow(res_grid))){
      x_block <- DelayedArray::read_block(x, x_grid[[coord2]])
      y_block <- DelayedArray::read_block(y, y_grid[[coord1]])
      res_block <- x_block %*% y_block
      DelayedArray::write_block(res_sink, res_grid[[coord2, coord1]], res_block)
    }
  }

  as(res_sink, "DelayedArray")
}


add_vector_to_each_column <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix + vector
}


add_vector_to_each_row <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if(length(vector) == 1){
    matrix + vector
  }else{
    t(t(matrix) + vector)
  }
}


subtract_vector_from_each_column <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix - vector
}


subtract_vector_from_each_row <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if(length(vector) == 1){
    matrix - vector
  }else{
    t(t(matrix) - vector)
  }
}


multiply_vector_to_each_column <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix * vector
}


multiply_vector_to_each_row <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if(length(vector) == 1){
    matrix * vector
  }else{
    t(t(matrix) * vector)
  }
}


divide_each_column_by_vector <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix / vector
}


divide_each_row_by_vector <- function(matrix, vector){
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if(length(vector) == 1){
    matrix / vector
  }else{
    t(t(matrix) / vector)
  }
}


