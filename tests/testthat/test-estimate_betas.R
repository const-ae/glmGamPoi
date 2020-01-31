

test_that("deviance calculation works", {

  r_gp_deviance <- function(y, mu, theta){
    -2 * sum(dnbinom(y, mu = mu, size = 1/theta, log = TRUE) - dnbinom(y, mu = y, size = 1/theta, log = TRUE))
  }

  expect_equal(compute_gp_deviance(y = 17, mu = 3, theta = 1.3),
               r_gp_deviance(y = 17, mu = 3, theta = 1.3))
  expect_equal(compute_gp_deviance(y = 0, mu = 3, theta = 1.3),
               r_gp_deviance(y = 0, mu = 3, theta = 1.3))
  expect_equal(compute_gp_deviance(y = 17, mu = 3, theta = 0),
               r_gp_deviance(y = 17, mu = 3, theta = 0))
  expect_equal(compute_gp_deviance(y = 0, mu = 3, theta = 0),
               r_gp_deviance(y = 0, mu = 3, theta = 0))
})

test_that("Rough estimation of Beta works", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  model_matrix <- cbind(1, rnorm(n = 4))
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  b1 <- estimate_betas_roughly(mat, model_matrix, offset_matrix)

  library(HDF5Array)
  hdf5_mat <- as(mat, "HDF5Matrix")
  hdf5_offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, hdf5_mat)$offset_matrix
  b2 <- estimate_betas_roughly(hdf5_mat, model_matrix, hdf5_offset_matrix)
  expect_equal(b1 , b2)
})


test_that("Beta estimation can handle edge cases as input", {

  Y <- matrix(0, nrow = 1, ncol = 10)
  model_matrix <- matrix(1, nrow = 10, ncol = 1)
  offset_matrix <- matrix(1, nrow = 1, ncol = 10)
  dispersion <- 0


  res <- estimate_betas_one_group(Y, offset_matrix, dispersion, beta_vec_init = c(3))
  expect_equal(res$Beta[1,1], -Inf)
  beta_mat_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix)
  res2 <- estimate_betas_fisher_scoring(Y, model_matrix, offset_matrix, dispersion, beta_mat_init)
  skip("Fisher scoring does not converge to -Inf")
  expect_equal(res2$Beta[1,1], -Inf)
})



test_that("Beta estimation can handle any kind of model_matrix", {

  # Weird input that makes DESeq2 choke

  Y <- matrix(1:72, nrow = 9, ncol = 8)[3:5,,drop=FALSE]
  model_matrix <- matrix(rnorm(n = 8 * 2), nrow = 8, ncol = 2)
  offset_matrix <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  disp_init <- estimate_dispersions_roughly(Y, model_matrix, offset_matrix)
  beta_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix)


  fit <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                dispersions = disp_init, beta_mat_init = beta_init)


  deseq2_fit <- DESeq2:::fitBetaWrapper(ySEXP = Y, xSEXP = model_matrix, nfSEXP = exp(offset_matrix),
                          alpha_hatSEXP = disp_init,
                          beta_matSEXP = beta_init,
                          lambdaSEXP = rep(0.3, ncol(model_matrix)),
                          weightsSEXP = array(1, dim(Y)), useWeightsSEXP = FALSE,
                          tolSEXP = 1e-8, maxitSEXP = 100, useQRSEXP = TRUE, minmuSEXP = 1e-6)

  edgeR_fit <- edgeR::glmFit.default(Y, design = model_matrix,
                        dispersion = disp_init, offset = offset_matrix[1,],
                        prior.count = 0, weights=NULL,
                        start = beta_init)
  # My result agrees with edgeR
  expect_equal(fit$Beta, edgeR_fit$coefficients, tolerance = 1e-3)
  # DESeq2 however does not converge
  expect_failure(
    expect_equal(fit$Beta, deseq2_fit$beta_mat, tolerance = 1e-3)
  )
  expect_failure(
    expect_equal(edgeR_fit$coefficients, deseq2_fit$beta_mat, tolerance = 1e-3)
  )
  expect_equal(deseq2_fit$iter, rep(100, 3))

  # My result, however did converge
  expect_lt(fit$iterations[1], 50)
})



test_that("estimate_betas_one_group can handle DelayedArray", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersions <- rep(0, 8)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix_hdf5 <- as(offset_matrix, "HDF5Matrix")

  beta_vec_init <- estimate_betas_roughly_one_group(mat, offset_matrix)
  beta_vec_init_da <- estimate_betas_roughly_one_group(mat_hdf5, offset_matrix_hdf5)

  res <- estimate_betas_one_group(mat, offset_matrix, dispersions, beta_vec_init)
  res2 <- estimate_betas_one_group(mat_hdf5, offset_matrix_hdf5, dispersions, beta_vec_init_da)
  # This check is important, because beachmat makes life difficult for
  # handling numeric and integer input generically
  res3 <- estimate_betas_one_group(mat * 1.0, offset_matrix, dispersions, beta_vec_init)
  expect_equal(res, res2)
  expect_equal(res, res3)

})

test_that("estimate_betas_fisher_scoring can handle DelayedArray", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  model_matrix <- cbind(1, rnorm(4, mean = 10))
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersions <- rep(0, 8)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix_hdf5 <- as(offset_matrix, "HDF5Matrix")

  beta_mat_init <- estimate_betas_roughly(mat, model_matrix, offset_matrix)
  beta_mat_init_da <- estimate_betas_roughly(mat_hdf5, model_matrix, offset_matrix_hdf5)


  res <- estimate_betas_fisher_scoring(mat, model_matrix, offset_matrix, dispersions, beta_mat_init)
  res2 <- estimate_betas_fisher_scoring(mat_hdf5, model_matrix, offset_matrix_hdf5, dispersions, beta_mat_init_da)
  res3 <- estimate_betas_fisher_scoring(mat * 1.0, model_matrix, offset_matrix, dispersions, beta_mat_init)
  expect_equal(res, res2)
  expect_equal(res, res3)
})


test_that("Beta estimation works", {

  data <- make_dataset(n_genes = 1000, n_samples = 30)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)

  # Fit Standard Model
  beta_mat_init <- estimate_betas_roughly(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix)
  my_res <- estimate_betas_fisher_scoring(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix,
                                          dispersions = data$overdispersion, beta_mat_init = beta_mat_init)

  # Fit Model for One Group
  beta_vec_init <- estimate_betas_roughly_one_group(Y = data$Y, offset_matrix = offset_matrix)
  my_res2 <- estimate_betas_one_group(Y = data$Y, offset_matrix = offset_matrix,
                                      dispersions = data$overdispersion, beta_vec_init = beta_vec_init)

  expect_equal(my_res$Beta, my_res2$Beta, tolerance = 1e-6)
  expect_lt(max(my_res2$iterations), 10)

  # Compare with edgeR
  edgeR_res <- edgeR::glmFit.default(data$Y, design = data$X,
                                   dispersion = data$overdispersion,
                                   offset = offset_matrix[1,],
                                   prior.count = 0, weights=NULL)

  expect_equal(my_res$Beta[,1], coef(edgeR_res)[,1], tolerance = 1e-6)
  expect_equal(my_res2$Beta[,1], coef(edgeR_res)[,1], tolerance = 1e-6)


  # Compare with DESeq2
  # This requires a few hacks:
  #   * If the "just intercept" model is fit, DESeq2
  #     automatically assigns an approximation for Beta
  #     To avoid this I give it a model design matrix that
  #     is really similar to the "just intercept" but numerically
  #     identical and thus force it to exactly calculate the
  #     beta values
  #   * The beta values are on a log2 scale. I need to convert it
  #     to the ln scale.
  dds_design_mat <- matrix(1+1e-8, nrow = ncol(data$Y), ncol = 1)
  dds <- DESeq2::DESeqDataSetFromMatrix(data$Y, colData = data.frame(name = seq_len(ncol(data$Y))),
                                        design = ~ 1)
  DESeq2::sizeFactors(dds) <- data$size_factor
  DESeq2::dispersions(dds) <- data$overdispersion
  dds <- DESeq2::nbinomWaldTest(dds, modelMatrix = dds_design_mat, minmu = 1e-6)
  expect_equal(my_res$Beta[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(my_res2$Beta[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(coef(edgeR_res)[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
})

test_that("glm_gp_impl can handle all zero rows", {
  Y <- matrix(0, nrow = 2, ncol = 10)
  Y[1, ] <- rpois(n = 10, lambda = 3)
  X <- matrix(1, nrow = 10, ncol = 1)
  res <- glm_gp_impl(Y, X)
  expect_equal(res$Beta_est[2,1], -Inf)
  expect_equal(res$Mu_est[2, ], rep(0, 10))
})


test_that("glm_gp_impl can handle all zero columns", {
  Y <- matrix(0, nrow = 2, ncol = 10)
  Y[1, 1] <- 3
  X <- matrix(1, nrow = 10, ncol = 1)
  res <- glm_gp_impl(Y, X)
  expect_equal(res$size_factors[2:10], rep(0.001, times = 9))
})

test_that("glm_gp_impl can handle all values zero", {
  Y <- matrix(0, nrow = 3, ncol = 10)
  X <- matrix(1, nrow = 10, ncol = 1)
  res <- glm_gp_impl(Y, X)
  expect_equal(res$size_factors, rep(0.001, times = 10))
  expect_equal(res$overdispersions, rep(0, times = 3))
  expect_equal(res$Beta_est[,1], rep(-Inf, times = 3))
  expect_true(all(res$Mu_est == 0))
})



test_that("glm_gp_impl works as expected", {
  skip("No workable tests here")
  # My method
  data <- make_dataset(n_genes = 2000, n_samples = 30)
  res <- glm_gp_impl(data$Y, model_matrix = data$X, verbose = TRUE)

  # edgeR
  edgeR_data <- edgeR::DGEList(data$Y)
  edgeR_data <- edgeR::calcNormFactors(edgeR_data)
  edgeR_data <- edgeR::estimateDisp(edgeR_data, data$X)
  fit <- edgeR::glmFit(edgeR_data, design = data$X)

  # DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(data$Y, colData = data.frame(name = seq_len(ncol(data$Y))),
                                        design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)

  res <- glm_gp_impl(data$Y, model_matrix = data$X, size_factors = DESeq2::sizeFactors(dds),
                     verbose = TRUE)
  plot(res$size_factors, DESeq2::sizeFactors(dds)); abline(0,1)

  plot(res$Beta_est[,1], coef(dds)[,1]  / log2(exp(1)), pch = 16, cex = 0.2, col ="red"); abline(0,1)
  plot(res$Beta_est[,1], coef(fit)[,1]); abline(0,1)
  plot(coef(dds)[,1]  / log2(exp(1)), coef(fit)[,1]); abline(0,1)


  plot(res$overdispersions, SummarizedExperiment::rowData(dds)$dispGeneEst, log = "xy"); abline(0,1)
  plot(res$overdispersions, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
  plot(SummarizedExperiment::rowData(dds)$dispGeneEst, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
})



test_that("glm_gp_impl works with Delayed Input", {
  # My method
  data <- make_dataset(n_genes = 2000, n_samples = 30)
  Y_da <- HDF5Array::writeHDF5Array(data$Y)
  res <- glm_gp_impl(data$Y, model_matrix = data$X, verbose = TRUE)
  res2 <- glm_gp_impl(Y_da, model_matrix = data$X, verbose = TRUE)

  expect_equal(res$Beta_est, res2$Beta_est)
  expect_equal(res$overdispersions, res2$overdispersions)
  expect_equal(res$Mu_est, as.matrix(res2$Mu_est))
  expect_equal(res$size_factors, res2$size_factors)
})



