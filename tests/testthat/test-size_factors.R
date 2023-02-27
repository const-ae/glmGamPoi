test_that("matrix vector operations work", {

  zmat <- matrix(0, nrow = 9, ncol = 4)
  onemat <- matrix(1, nrow = 9, ncol = 4)
  expect_equal(add_vector_to_each_column(zmat, 1:9)[,1],
               1:9)
  expect_equal(subtract_vector_from_each_column(zmat, 1:9)[,1],
               -(1:9))
  expect_equal(multiply_vector_to_each_column(onemat, 1:9)[,1],
               1:9)
  expect_equal(divide_each_column_by_vector(onemat, 1:9)[,1],
               1/(1:9))

  expect_equal(add_vector_to_each_row(zmat, 1:4)[1, ],
               1:4)
  expect_equal(subtract_vector_from_each_row(zmat, 1:4)[1, ],
               -(1:4))
  expect_equal(multiply_vector_to_each_row(onemat, 1:4)[1, ],
               1:4)
  expect_equal(divide_each_row_by_vector(onemat, 1:4)[1, ],
               1/(1:4))

})



test_that("estimate size factor works", {

  mat <- matrix(rpois(n = 32, lambda = 5), nrow = 8, ncol = 4)
  library(HDF5Array)
  hdf5_mat <- as(mat, "HDF5Matrix")

  expect_equal(estimate_size_factors(mat, method = "deconvolution"),
               estimate_size_factors(hdf5_mat, method = "deconvolution"))

  expect_equal(estimate_size_factors(mat, method = "normed_sum"),
               estimate_size_factors(hdf5_mat, method = "normed_sum"))

  expect_equal(estimate_size_factors(mat, method = "poscounts"),
               estimate_size_factors(hdf5_mat, method = "poscounts"))

  expect_equal(estimate_size_factors(mat, method = "ratio"),
               estimate_size_factors(hdf5_mat, method = "ratio"))

  sf <- estimate_size_factors(mat, method = "ratio")
  sf_deseq <- DESeq2::estimateSizeFactorsForMatrix(mat)
  sf_deseq <- sf_deseq / exp(mean(log(sf_deseq)))
  expect_equal(sf, sf_deseq)

  mat <- matrix(0, nrow = 8, ncol = 4)
  mat[1,1] <- 7
  mat[3, 1] <- 7
  hdf5_mat <- as(mat, "HDF5Matrix")

  # Throws an error
  # expect_equal(estimate_size_factors(mat, method = "deconvolution"),
  #              estimate_size_factors(hdf5_mat, method = "deconvolution"))

  expect_equal(estimate_size_factors(mat, method = "normed_sum"),
               estimate_size_factors(hdf5_mat, method = "normed_sum"))

  expect_equal(estimate_size_factors(mat, method = "poscounts"),
               estimate_size_factors(hdf5_mat, method = "poscounts"))

  expect_error(estimate_size_factors(mat, method = "ratio"))
  expect_error(estimate_size_factors(hdf5_mat, method = "asdf"))

})



test_that("combine_size_factors_and_offset works", {
  mat <- matrix(1:32, nrow = 8, ncol = 4)
  offset_num <- 0
  sf <- 1.1

  off_and_sf <- combine_size_factors_and_offset(offset_num, sf, mat)
  expect_equal(off_and_sf$size_factors, rep(1.1, times = 4))
  expect_equal(off_and_sf$offset_matrix[1,1], log(1.1))
  expect_is(off_and_sf$offset_matrix, "matrix")

  sf2 <- exp(rnorm(n = 4))

  off_and_sf <- combine_size_factors_and_offset(offset_num, sf2, mat)
  expect_equal(off_and_sf$size_factors, sf2)
  expect_equal(off_and_sf$offset_matrix[1,], log(sf2))
  expect_is(off_and_sf$offset_matrix, "matrix")


  sf3 <- TRUE
  est_sf <- estimate_size_factors(mat, method = "deconvolution")
  off_and_sf <- combine_size_factors_and_offset(offset_num, sf3, mat)
  expect_equal(off_and_sf$size_factors, est_sf)
  expect_equal(off_and_sf$offset_matrix[1,], log(est_sf))
  expect_is(off_and_sf$offset_matrix, "matrix")

  library(HDF5Array)
  hdf5_mat <- as(mat, "HDF5Matrix")
  off_and_sf <- combine_size_factors_and_offset(offset_num, sf3, hdf5_mat)
  expect_equal(off_and_sf$size_factors, est_sf)
  expect_equal(off_and_sf$offset_matrix[1,], log(est_sf))
  expect_is(off_and_sf$offset_matrix, "HDF5Matrix")

  offset_num2 <- 1:4
  off_and_sf <- combine_size_factors_and_offset(offset_num2, est_sf, hdf5_mat)
  expect_equal(off_and_sf$size_factors, est_sf)
  expect_equal(off_and_sf$offset_matrix[1,], 1:4 + log(est_sf))
  expect_is(off_and_sf$offset_matrix, "HDF5Matrix")
})


