

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
  dispersion <- rep(0, 10)


  res <- estimate_betas_one_group(Y, offset_matrix, dispersion, beta_vec_init = c(3))
  expect_equal(res$Beta[1,1], -Inf)
  skip("Fisher scoring does not converge to -Inf")
  res2 <- estimate_betas(Y, model_matrix, offset_matrix, dispersion)
  expect_equal(res2$Beta[1,1], -Inf)
})


test_that("estimate_betas_one_group can handle DelayedArray", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersion <- rep(0, 10)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix_hdf5 <- as(offset_matrix, "HDF5Matrix")

  res <- estimate_betas_one_group(mat, offset_matrix, dispersion)
  res2 <- estimate_betas_one_group(mat_hdf5, offset_matrix_hdf5, dispersion)
  # This check is important, because beachmat makes life difficult for
  # handling numeric and integer input generically
  res3 <- estimate_betas_one_group(mat * 1.0, offset_matrix, dispersion)
  expect_equal(res, res2)
  expect_equal(res, res3)

})

test_that("estimate_betas_fisher_scoring can handle DelayedArray", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  design_matrix <- cbind(1, rnorm(4, mean = 10))
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersion <- rep(0, 10)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix_hdf5 <- as(offset_matrix, "HDF5Matrix")

  res <- estimate_betas_fisher_scoring(mat, design_matrix, offset_matrix, dispersion)
  res2 <- estimate_betas_fisher_scoring(mat_hdf5, design_matrix, offset_matrix_hdf5, dispersion)
  res3 <- estimate_betas_fisher_scoring(mat * 1.0, design_matrix, offset_matrix, dispersion)
  expect_equal(res, res2)
  expect_equal(res, res3)
})


test_that("Beta estimation works", {

  data <- make_dataset(n_genes = 1000, n_samples = 30)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)

  # Fit Standard Model
  my_res <- estimate_betas_fisher_scoring(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix,
                                          dispersion = data$overdispersion)

  # Fit Model for One Group
  my_res2 <- estimate_betas_one_group(Y = data$Y, offset_matrix = offset_matrix,
                                      dispersion = data$overdispersion)

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
  Y <- matrix(0, nrow = 1, ncol = 10)
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
  res <- glm_gp_impl(data$Y, design_matrix = data$X, verbose = TRUE)

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

  res <- glm_gp_impl(data$Y, design_matrix = data$X, size_factors = DESeq2::sizeFactors(dds),
                     verbose = TRUE)
  plot(res$size_factors, DESeq2::sizeFactors(dds)); abline(0,1)

  plot(res$Beta_est[,1], coef(dds)[,1]  / log2(exp(1)), pch = 16, cex = 0.2, col ="red"); abline(0,1)
  plot(res$Beta_est[,1], coef(fit)[,1]); abline(0,1)
  plot(coef(dds)[,1]  / log2(exp(1)), coef(fit)[,1]); abline(0,1)


  plot(res$overdispersions, SummarizedExperiment::rowData(dds)$dispGeneEst, log = "xy"); abline(0,1)
  plot(res$overdispersions, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
  plot(SummarizedExperiment::rowData(dds)$dispGeneEst, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
})







## More Benchmark:
# data <- make_dataset(n_genes = 1000, n_samples = 3000)
# profvis::profvis(
#   glm_gp_impl(data$Y, design_matrix = data$X, verbose = TRUE)
# )
# bench::mark(
#   glmGamPoi = {
#     glm_gp_impl(data$Y, design_matrix = data$X, verbose = TRUE)
#   # }, DESeq2 = {
#   #   dds <- DESeq2::DESeqDataSetFromMatrix(data$Y, colData = data.frame(name = seq_len(ncol(data$Y))),
#   #                                         design = ~ 1)
#   #   dds <- DESeq2::estimateSizeFactors(dds)
#   #   dds <- DESeq2::estimateDispersions(dds)
#   #   dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
#   }, edgeR = {
#     edgeR_data <- edgeR::DGEList(data$Y)
#     edgeR_data <- edgeR::calcNormFactors(edgeR_data)
#     edgeR_data <- edgeR::estimateDisp(edgeR_data, data$X)
#     fit <- edgeR::glmFit(edgeR_data, design = data$X)
#   }, check = FALSE
# )


## Benchmark:
# bench::mark(
#   gp_general = estimate_betas(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix,
#                            dispersions = data$overdispersion),
#   gp_one_group = estimate_betas_one_group(Y = data$Y, offset_matrix = offset_matrix,
#                                       dispersions = data$overdispersion),
#   edgeR_one_group = edgeR::glmFit.default(data$Y, design = data$X,
#                                 dispersion = data$overdispersion,
#                                 offset = offset_matrix[1,],
#                                 prior.count = 0, weights=NULL),
#   edgeR_general = edgeR::glmFit.default(data$Y,
#                                 design = matrix(1+rnorm(n=ncol(data$Y), mean=0, sd=1e-8),
#                                                 nrow = ncol(data$Y), ncol = 1),
#                                 dispersion = data$overdispersion,
#                                 offset = offset_matrix[1,],
#                                 prior.count = 0, weights=NULL),
#   DESeq2_one_group = DESeq2:::fitNbinomGLMs(dds, modelFormula = ~ 1, alpha_hat = data$overdispersion,
#                                          lambda = 1e-6, renameCols = FALSE, minmu = 1e-6),
#   DESeq2_general =  DESeq2:::fitNbinomGLMs(dds, modelMatrix = matrix(1+1e-8, nrow = ncol(data$Y), ncol = 1),
#                                          alpha_hat = data$overdispersion,
#                                          lambda = 1e-6, renameCols = FALSE, minmu = 1e-6),
#   check = FALSE
# ) %>%
#   arrange(median) %>%
#   mutate(relative = c(median / first(median)))







