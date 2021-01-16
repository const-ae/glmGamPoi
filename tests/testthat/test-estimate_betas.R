

setup({
  # Limit the number of workers so CRAN is happy
  DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers = 2))
})

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


  res <- estimate_betas_group_wise(Y, offset_matrix, dispersion, beta_group_init = matrix(3, nrow = 1, ncol = 1), groups = 1, model_matrix = model_matrix)
  expect_equal(res$Beta[1,1], -1e8)
  beta_mat_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix)
  res2 <- estimate_betas_fisher_scoring(Y, model_matrix, offset_matrix, dispersion, beta_mat_init)
  expect_lt(res2$Beta[1,1], -15)
})

test_that("Beta estimation can handle edge case (2)", {
  # In version 1.1.13, glmGamPoi converged to an extreme result (mu = 1e50) for this input.
  # With the introduction of the max_rel_mu_change parameter, this seems to be fixed

  y <- c(0, 0, 14, 2, 0, 0, 0, 0, 10, 12, 6, 2, 0, 4, 1, 1, 2, 6, 165, 2, 1, 0, 0, 0, 259, 2050, 715, 0, 0, 96, 2658, 149, 56, 7, 0, 0, 0, 0, 0, 0, 0, 5, 9, 1, 1, 0, 0, 1, 1, 5, 7, 0, 0, 1, 0, 3, 6, 19, 29, 0, 0, 0, 1, 4, 73, 3, 4, 1, 0, 1, 0, 0, 5, 169, 58, 1, 0, 32, 0, 2, 1, 1, 170, 30, 0, 1, 0, 4, 123, 1655, 1292, 101, 0, 732, 2866, 207, 3, 6, 3, 0, 0, 2, 0, 1, 0, 110, 27, 0, 0, 0, 0, 1, 51, 3, 198, 1, 0, 1, 0, 78, 9, 2, 142, 2, 0, 2, 1, 1)
  clust <- c("B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells",
             "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cells", "B cells", "CD14+ Monocytes", "CD4 T cells", "CD8 T cells", "Dendritic cells", "FCGR3A+ Monocytes", "Megakaryocytes", "NK cell")
  cond <- c("ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim", "stim")
  sf <- c(1.25495746, 2.84060391, 3.54532248, 0.74403676, 0.12168579, 1.29108360, 0.12000349, 0.65455177, 5.44150278, 12.27749619, 10.45289957, 1.95654869, 0.13524118, 3.82396827, 0.31255044, 1.73958695, 1.48917348, 6.17830352, 4.72301313, 5.18268863, 0.09076828, 1.82555638, 0.08546259, 1.12790128, 0.39478863, 1.69810120, 2.50113680, 0.26218952, 0.05861062, 0.29934013, 0.36996404, 0.23265235, 0.55359268, 3.64881579, 2.02365705, 0.20378379, 0.02849112, 0.38153519, 0.08851444, 0.37638048, 1.51409512, 6.09265761, 13.59420239, 0.64251691, 0.38507232, 0.67384420, 0.20171328, 1.01657884, 2.55432311, 5.04959345, 12.65532171, 1.37339728, 0.16297743, 0.65830457, 0.20081821, 2.34212786, 2.65519592, 4.60271847, 14.42841431, 0.61404735, 0.33373006, 1.61952954, 0.20445239, 1.01918855, 1.71367318, 2.73565477, 4.88628172, 1.19830951, 0.23995307, 1.79371146, 0.10530501, 1.02351290, 4.40897451, 9.49905081, 9.37883164, 1.63402313, 0.17510934, 2.88147498, 0.20539059, 2.25509082, 1.48036301, 5.38691456, 4.21155324, 5.01258303, 0.19591153, 2.02167280, 0.11638010, 1.98069390, 0.57433016, 2.00913110, 3.63831225, 0.33962887, 0.20308284, 0.39905907, 0.46129308, 0.29367857, 0.69204747, 2.43350004, 2.24865281, 0.16404503, 0.11647715, 0.52166148, 0.04277982, 0.48455401, 1.23137302, 4.23807090, 11.94627878, 0.39999727, 0.23969425, 0.29627750, 0.14492514, 1.20180350, 2.62172262, 4.84203529, 12.69336739, 1.23926685, 0.28333679, 1.02518441, 0.24187261, 2.28643968, 3.42665620, 5.63192529, 19.61644098, 0.48511477, 0.52967394, 2.41232041, 0.37996074, 1.65937613)

  model_matrix <- model.matrix(~ cond + clust)
  offset_matrix <- add_vector_to_each_row(matrix(0, nrow = 1, ncol = length(y)), log(sf))
  Y <- matrix(y, nrow = 1)
  disp_init <- estimate_dispersions_roughly(Y, model_matrix, offset_matrix = offset_matrix)
  beta_init <- estimate_betas_roughly(Y, model_matrix, offset_matrix = offset_matrix)
  beta_res <- estimate_betas_fisher_scoring(Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                                        dispersions = disp_init, beta_mat_init = beta_init)
  beta_res
  expect_lt(beta_res$deviances, 100)
  expect_true(all(calculate_mu(beta_res$Beta, model_matrix, offset_matrix) < 1e5))
  # betaRes <- fitBeta_fisher_scoring(Y, model_matrix, exp(offset_matrix), thetas = disp_init, beta_matSEXP = beta_init,
  #                                   ridge_penalty = 1e-6, tolerance = 1e-8, max_iter =  1000)
  # betaRes
})


test_that("Groupwise beta estimation works", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  b1 <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = 1)
  b2 <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = rep(1, times = ncol(mat)))
  expect_equal(b1, b2)

  model_matrix <- cbind(c(1,1,0,0), c(0,0,1,1))
  b3 <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = c(1, 1, 2, 2))
  b4 <- cbind(
    estimate_betas_roughly_group_wise(mat[,1:2], offset_matrix[,1:2], groups = 1),
    estimate_betas_roughly_group_wise(mat[,3:4], offset_matrix[,3:4], groups = 1)
  )
  expect_equal(b3, b4)


  b5 <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = c(1, 1, 2, 2))
  b6 <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = c(2, 2, 1, 1))
  expect_equal(b5, b6)
})




test_that("Beta estimation can handle any kind of model_matrix", {
  skip_if_not_installed("DESeq2")
  # Weird input that makes DESeq2 choke
  set.seed(1)
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

test_that("estimate_betas_group_wise properly rescales result", {

  dat <- make_dataset(n_genes = 20, n_samples = 30)
  mat <- dat$Y
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersions <- dat$overdispersion

  df <- data.frame(cat1 = sample(LETTERS[1:3], size = 30, replace = TRUE))
  model_matrix <- model.matrix(~ cat1, data = df)

  groups <- get_groups_for_model_matrix(model_matrix)
  beta_group_init <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = groups)

  res <- estimate_betas_group_wise(mat, offset_matrix, dispersions, beta_group_init, groups = groups, model_matrix = model_matrix)

  model_matrix2 <- model_matrix * 3
  res2 <- estimate_betas_group_wise(mat, offset_matrix, dispersions, beta_group_init = beta_group_init, groups = groups, model_matrix = model_matrix2)
  expect_equal(res$Beta, res2$Beta * 3)
  expect_equal(calculate_mu(res$Beta, model_matrix, offset_matrix),
               calculate_mu(res2$Beta, model_matrix2, offset_matrix))
})

test_that("estimate_betas_group_wise can handle extreme case", {

  y <- matrix(c(1, rep(0, 500)), nrow = 1)
  res <- fitBeta_one_group(y, offset_matrix = matrix(0, nrow = 1, ncol = 501), thetas = 2.1, beta_start_values = -10, tolerance = 1e-8,  maxIter = 100)
  expect_false(res$iter == 100)
  expect_false(is.na(res$beta))

})


test_that("estimate_betas_group_wise can handle DelayedArray", {

  mat <- matrix(1:32, nrow = 8, ncol = 4)
  offset_matrix <- combine_size_factors_and_offset(0, size_factors = TRUE, mat)$offset_matrix
  dispersions <- rep(0, 8)
  model_matrix <- matrix(1, nrow = 4)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix_hdf5 <- as(offset_matrix, "HDF5Matrix")


  beta_vec_init <- estimate_betas_roughly_group_wise(mat, offset_matrix, groups = 1)
  beta_vec_init_da <- estimate_betas_roughly_group_wise(mat_hdf5, offset_matrix_hdf5, groups = 1)

  res <- estimate_betas_group_wise(mat, offset_matrix, dispersions, beta_vec_init, groups = 1, model_matrix = model_matrix)
  res2 <- estimate_betas_group_wise(mat_hdf5, offset_matrix_hdf5, dispersions, beta_vec_init_da, groups = 1, model_matrix = model_matrix)
  # This check is important, because beachmat makes life difficult for
  # handling numeric and integer input generically
  res3 <- estimate_betas_group_wise(mat * 1.0, offset_matrix, dispersions, beta_vec_init, groups = 1, model_matrix = model_matrix)
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
  skip_if_not(is_macos(), "Beta estimation is unprecise on Non-MacOS architectures")
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")
  data <- make_dataset(n_genes = 1000, n_samples = 30)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)
  model_matrix <- matrix(1.1, nrow = 30)

  # Fit Standard Model
  beta_mat_init <- estimate_betas_roughly(Y = data$Y, model_matrix = model_matrix, offset_matrix = offset_matrix)
  my_res <- estimate_betas_fisher_scoring(Y = data$Y, model_matrix = model_matrix, offset_matrix = offset_matrix,
                                          dispersions = data$overdispersion, beta_mat_init = beta_mat_init)

  # Fit Model for One Group
  beta_vec_init <- estimate_betas_roughly_group_wise(Y = data$Y, offset_matrix = offset_matrix, groups = 1)
  my_res2 <- estimate_betas_group_wise(Y = data$Y, offset_matrix = offset_matrix,
                                      dispersions = data$overdispersion, beta_group_init = beta_vec_init, groups = 1, model_matrix = model_matrix)

  expect_equal(my_res$Beta, my_res2$Beta, tolerance = 1e-6)
  expect_lt(max(my_res2$iterations), 10)

  # Compare with edgeR
  edgeR_res <- edgeR::glmFit.default(data$Y, design = model_matrix,
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
  dds <- DESeq2::DESeqDataSetFromMatrix(data$Y, colData = data.frame(name = seq_len(ncol(data$Y))),
                                        design = ~ 1)
  DESeq2::sizeFactors(dds) <- data$size_factor
  DESeq2::dispersions(dds) <- data$overdispersion
  dds <- DESeq2::nbinomWaldTest(dds, modelMatrix = model_matrix, minmu = 1e-6)
  expect_equal(my_res$Beta[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(my_res2$Beta[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(coef(edgeR_res)[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
})



test_that("Fisher scoring and diagonal fisher scoring give consistent results", {
  set.seed(1)
  data <- make_dataset(n_genes = 1, n_samples = 3000)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)

  # Fit Standard Model
  beta_mat_init <- estimate_betas_roughly(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix)
  res1 <- fitBeta_fisher_scoring(Y = data$Y, model_matrix = data$X, exp_offset_matrix = exp(offset_matrix),
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = 0,
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  res2 <- fitBeta_diagonal_fisher_scoring(Y = data$Y, model_matrix = data$X, exp_offset_matrix = exp(offset_matrix),
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init,
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  expect_equal(res1, res2  )

  set.seed(1)
  df <- data.frame(
    city = sample(c("Heidelberg", "Berlin", "New York"), size = 3000, replace = TRUE),
    fruit = sample(c("Apple", "Cherry", "Banana"), size = 3000, replace = TRUE),
    age = rnorm(3000, mean = 50, sd = 15),
    car = sample(c("blue", "big", "truck"), size = 3000, replace = TRUE),
    letter = LETTERS[1:10]
  )
  new_model_matrix <- model.matrix(~ . - 1, df)
  beta_mat_init <- estimate_betas_roughly(Y = data$Y, model_matrix = new_model_matrix, offset_matrix = offset_matrix)
  res1 <- fitBeta_fisher_scoring(Y = data$Y, model_matrix = new_model_matrix, exp_offset_matrix = exp(offset_matrix),
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = rep(0, ncol(new_model_matrix)),
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  res2 <- fitBeta_diagonal_fisher_scoring(Y = data$Y, model_matrix = new_model_matrix, exp_offset_matrix = exp(offset_matrix),
                                          thetas = data$overdispersion, beta_matSEXP = beta_mat_init,
                                          tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  1000)
  expect_gt(cor(c(res1$beta_mat), c(res2$beta_mat)), 0.99)
  expect_gt(res2$iter, res1$iter)
})


test_that("Fisher scoring and ridge penalized fisher scoring give consistent results", {

  data <- make_dataset(n_genes = 1, n_samples = 3000)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)

  # Fit Standard Model
  beta_mat_init <- estimate_betas_roughly(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix)
  res1 <- fitBeta_fisher_scoring(Y = data$Y, model_matrix = data$X, exp_offset_matrix = exp(offset_matrix),
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = 0,
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  res2 <- fitBeta_fisher_scoring(Y = data$Y, model_matrix = data$X, exp_offset_matrix = exp(offset_matrix),
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = 1e-10,
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  expect_equal(res1, res2)

  set.seed(1)
  size <- 25
  df <- data.frame(
    city = sample(c("Heidelberg", "Berlin", "New York"), size = size, replace = TRUE),
    fruit = sample(c("Apple", "Cherry", "Banana"), size = size, replace = TRUE),
    age = rnorm(size, mean = 50, sd = 1e-5),
    car = sample(c("blue", "big", "truck"), size = size, replace = TRUE)
  )
  new_model_matrix <- model.matrix(~ . - 1, df)
  beta_mat_init <- estimate_betas_roughly(Y = data$Y[,1:size,drop=FALSE], model_matrix = new_model_matrix, offset_matrix = offset_matrix[,1:size,drop=FALSE])
  res1 <- fitBeta_fisher_scoring(Y = data$Y[,1:size,drop=FALSE], model_matrix = new_model_matrix, exp_offset_matrix = exp(offset_matrix)[,1:size,drop=FALSE],
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = rep(0, ncol(new_model_matrix)),
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  res2 <- fitBeta_fisher_scoring(Y = data$Y[,1:size,drop=FALSE], model_matrix = new_model_matrix, exp_offset_matrix = exp(offset_matrix)[,1:size,drop=FALSE],
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = rep(1e-30, ncol(new_model_matrix)),
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  res3 <- fitBeta_fisher_scoring(Y = data$Y[,1:size,drop=FALSE], model_matrix = new_model_matrix, exp_offset_matrix = exp(offset_matrix)[,1:size,drop=FALSE],
                                 thetas = data$overdispersion, beta_matSEXP = beta_mat_init, ridge_penalty = rep(50, ncol(new_model_matrix)),
                                 tolerance = 1e-8, max_rel_mu_change = 1e5, max_iter =  100)
  expect_equal(res1, res2, tolerance = 1e-6)
  expect_lt(res3$beta_mat[6], res1$beta_mat[6])  # The age column is much smaller
})


test_that("glm_gp_impl can handle all zero rows", {
  Y <- matrix(0, nrow = 2, ncol = 10)
  Y[1, ] <- rpois(n = 10, lambda = 3)
  X <- matrix(1, nrow = 10, ncol = 1)
  res <- glm_gp_impl(Y, X)
  expect_equal(res$Beta[2,1], -1e8)
  expect_equal(res$Mu[2, ], rep(0, 10))
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
  expect_equal(res$Beta[,1], rep(-1e8, times = 3))
  expect_true(all(res$Mu == 0))
})


test_that("glm_gp_impl can handle dispersion of zero", {
  Y <- matrix(rnbinom(10, mu = 10, size = 1 / 2.3), nrow = 1, ncol = 10)
  X <- cbind(1, rnorm(10))
  res <- glm_gp_impl(Y, X, overdispersion = 0, size_factors = FALSE)
  res2 <- glm(c(Y) ~ X - 1, family = "poisson")

  expect_equal(c(res$Beta), unname(coef(res2)))

})


test_that("glm_gp_impl can handle weird input", {
  Y <- matrix(c(7, 0, 0, 0, 0), nrow = 1)
  X <- cbind(1, c(0.64, -2.7, -0.94, 0.6, 0.56))
  offset <- matrix(0, nrow = 1, ncol = 5)
  init <- matrix(c(1,1), nrow = 1)
  # This used to return c(NA, NA) because mu got exactly zero
  res <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init)
  # fitBeta_diagonal_fisher_scoring(Y, X, exp(offset), 0, init, tolerance = 1e-8, max_iter =  5000000)
  expect_false(any(is.na(c(res$Beta))))
})


test_that("glm_gp_impl can handle weird input 2", {
  Y <- matrix(c(34, 130, 1, 27, 1), nrow = 1)
  X <- cbind(c(0.03, -0.4, -0.4, 0.8, 0.1),
             c(-0.1, -0.7, 0.7, -0.03, 0.2))
  offset <- matrix(0, nrow = 1, ncol = 5)
  init <- matrix(c(3000, -141), nrow = 1)
  res <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init)
  expect_false(any(is.na(c(res$Beta))))
})




# test_that("glm_gp_impl works as expected", {
#   skip("No workable tests here")
#   # My method
#   data <- make_dataset(n_genes = 2000, n_samples = 30)
#   res <- glm_gp_impl(data$Y, model_matrix = data$X, verbose = TRUE)
#
#   # edgeR
#   edgeR_data <- edgeR::DGEList(data$Y)
#   edgeR_data <- edgeR::calcNormFactors(edgeR_data)
#   edgeR_data <- edgeR::estimateDisp(edgeR_data, data$X)
#   fit <- edgeR::glmFit(edgeR_data, design = data$X)
#
#   # DESeq2
#   dds <- DESeq2::DESeqDataSetFromMatrix(data$Y, colData = data.frame(name = seq_len(ncol(data$Y))),
#                                         design = ~ 1)
#   dds <- DESeq2::estimateSizeFactors(dds)
#   dds <- DESeq2::estimateDispersions(dds)
#   dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
#
#   res <- glm_gp_impl(data$Y, model_matrix = data$X, size_factors = DESeq2::sizeFactors(dds),
#                      verbose = TRUE)
#   plot(res$size_factors, DESeq2::sizeFactors(dds)); abline(0,1)
#
#   plot(res$Beta[,1], coef(dds)[,1]  / log2(exp(1)), pch = 16, cex = 0.2, col ="red"); abline(0,1)
#   plot(res$Beta[,1], coef(fit)[,1]); abline(0,1)
#   plot(coef(dds)[,1]  / log2(exp(1)), coef(fit)[,1]); abline(0,1)
#
#
#   plot(res$overdispersions, SummarizedExperiment::rowData(dds)$dispGeneEst, log = "xy"); abline(0,1)
#   plot(res$overdispersions, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
#   plot(SummarizedExperiment::rowData(dds)$dispGeneEst, edgeR_data$tagwise.dispersion, log = "xy"); abline(0,1)
# })



test_that("glm_gp_impl works with Delayed Input", {
  # My method
  data <- make_dataset(n_genes = 2000, n_samples = 30)
  Y_da <- HDF5Array::writeHDF5Array(data$Y)
  res <- glm_gp_impl(data$Y, model_matrix = data$X, verbose = TRUE)
  res2 <- glm_gp_impl(Y_da, model_matrix = data$X, verbose = TRUE)

  expect_equal(res$Beta, res2$Beta)
  expect_equal(res$overdispersions, res2$overdispersions)
  expect_equal(res$Mu, as.matrix(res2$Mu))
  expect_equal(res$size_factors, res2$size_factors)
})



test_that("ridge penalization works as expected", {
  set.seed(1)
  X <- cbind(1, rnorm(n = 100))
  Y <- matrix(rpois(n = nrow(X), lambda = 4), nrow = 1)
  offset <- matrix(0, nrow = 1, ncol = nrow(X))
  init <- matrix(c(1,1), nrow = 1)
  # This used to return c(NA, NA) because mu got exactly zero
  res1 <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init,
                                       ridge_penalty = 0)
  res2 <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init,
                                       ridge_penalty = c(0, 10))
  res1
  res2

  expect_lt(res2$Beta[2], res1$Beta[2])

  # Group estimator
  X <- cbind(rep(c(0, 1, 1), length = 100), rep(c(1, 0, 0), length = 100))
  glm_gp(Y ~ X - 1, ridge_penalty = 10, verbose = TRUE)

  res1 <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init,
                                        ridge_penalty = NULL)
  res2 <- estimate_betas_fisher_scoring(Y, X, offset, dispersions = 0, beta_mat_init = init,
                                        ridge_penalty = c(1, 1))

  res3 <- estimate_betas_group_wise(Y, offset, dispersions = 0, beta_mat_init = init,
                                    groups = rep(c(1,2,2), length = 100), model_matrix = X)


  expect_equal(res1[c(1,3)], res3[c(1,3)])
  expect_gt(res2$deviances, res1$deviances)
  expect_lt(res2$Beta[1], res1$Beta[1])

})




