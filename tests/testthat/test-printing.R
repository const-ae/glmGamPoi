skip("Printing tests should be run interactively")

test_that("summary of matrix looks good", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  cat(format_matrix(mat))
  colnames(mat) <- LETTERS[1:4]
  cat(format_matrix(mat))

  mat <- matrix(1:4, nrow = 1, ncol = 4)
  cat(format_matrix(mat))
  colnames(mat) <- LETTERS[1:4]
  cat(format_matrix(mat))


  mat <- matrix(numeric(0), nrow = 0, ncol = 4)
  cat(format_matrix(mat))
  colnames(mat) <- LETTERS[1:4]
  cat(format_matrix(mat))

  mat <- matrix(numeric(0), nrow = 4, ncol = 0)
  cat(format_matrix(mat))
  rownames(mat) <- LETTERS[1:4]
  cat(format_matrix(mat))


})


test_that("summary of matrix can handle NAs", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  mat[1, 1] <- NA
  cat(format_matrix(my_matrix_summary(mat)))
  colnames(mat) <- LETTERS[1:4]
  cat(format_matrix(my_matrix_summary(mat)))

})


test_that("fit printing works", {
  mat <- matrix(1:32, nrow = 8, ncol = 4)
  fit <- glm_gp(mat)
  print(fit)
  print(summary(fit))

  fit <- glm_gp(mat[1:3, ,drop=FALSE], size_factors = FALSE)
  print(fit)
  print(summary(fit))

  model_matrix <- matrix(rnorm(n = 32), nrow = 8, ncol = 4)
  mat <- matrix(1:72, nrow = 9, ncol = 8)
  fit <- glm_gp(mat, design = model_matrix)
  print(fit)
  print(summary(fit))

  fit <- glm_gp(mat[1:3, ], design = model_matrix[, 1:2])
  print(fit)
  print(summary(fit))

})


test_that("fit printing can handle DelayedArrays", {
  mat <- matrix(1:32, nrow = 8, ncol = 4)
  mat5 <- HDF5Array::writeHDF5Array(mat)
  fit <- glm_gp(mat5)
  print(fit)
  print(summary(fit))
})


