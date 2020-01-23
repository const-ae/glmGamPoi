

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


