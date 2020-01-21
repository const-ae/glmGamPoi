

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


test_that("Beta estimation works", {

  data <- make_dataset(n_genes = 1000, n_samples = 30)
  offset_matrix <- matrix(log(data$size_factor), nrow=nrow(data$Y), ncol = ncol(data$Y), byrow = TRUE)

  # Fit Standard Model
  my_res <- estimate_betas(Y = data$Y, model_matrix = data$X, offset_matrix = offset_matrix,
                          dispersions = data$overdispersion)

  # Fit Model for One Group
  my_res2 <- estimate_betas_one_group(Y = data$Y, offset_matrix = offset_matrix,
                                      dispersions = data$overdispersion)

  expect_equal(my_res$beta_mat[,1], my_res2$beta, tolerance = 1e-6)
  expect_lt(max(my_res2$iter), 10)

  # Compare with edgeR
  edgeR_res <- edgeR::glmFit.default(data$Y, design = data$X,
                                   dispersion = data$overdispersion,
                                   offset = offset_matrix[1,],
                                   prior.count = 0, weights=NULL)

  expect_equal(my_res$beta_mat[,1], coef(edgeR_res)[,1], tolerance = 1e-6)
  expect_equal(my_res2$beta, coef(edgeR_res)[,1], tolerance = 1e-6)


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
  expect_equal(my_res$beta_mat[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(my_res2$beta, coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
  expect_equal(coef(edgeR_res)[,1], coef(dds)[,1] / log2(exp(1)), tolerance = 1e-6)
})



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







