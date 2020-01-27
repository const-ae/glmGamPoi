


test_that("DelayedMultiplication works properly", {

  Beta <- matrix(c(rnorm(n = 200), rnorm(n = 200, mean = 10)), nrow = 200, ncol = 2)
  model_matrix <- cbind(1, rpois(n = 7000, lambda = 3))

  res <- Beta %*% t(model_matrix)

  Beta_da <- DelayedArray::DelayedArray(Beta)
  model_matrix_da <- DelayedArray::DelayedArray(model_matrix)
  res2 <- delayed_matrix_multiply(Beta_da, t(model_matrix_da))

  expect_equal(res, as.matrix(res2))


})


test_that("calculating mu with DelayedArrays works properly", {

  Beta <- matrix(c(rnorm(n = 20), rnorm(n = 20, mean = 10)), nrow = 20, ncol = 2)
  model_matrix <- cbind(1, rpois(n = 300, lambda = 3))
  offset_matrix <- matrix(1.3, nrow = 20, ncol = 300)
  mu <- calculate_mu(Beta, model_matrix, offset_matrix)

  Beta_da <- DelayedArray::DelayedArray(Beta)
  model_matrix_da <- DelayedArray::DelayedArray(model_matrix)
  offset_matrix_hdf5 <- HDF5Array::writeHDF5Array(offset_matrix)
  mu2 <- calculate_mu(Beta_da, model_matrix_da, offset_matrix_hdf5)

  expect_equal(mu, as.matrix(mu2))
})


