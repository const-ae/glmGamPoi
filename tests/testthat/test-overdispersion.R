


test_that("makeCumSumLookupVector works", {
  vector <- rpois(n = 100, lambda = 10)
  cslv <- makeCumSumLookupVector(vector)
  lv <- table(factor(vector, levels = seq_len(max(vector)+1)-1, ordered = TRUE))
  cslv2 <- rev(cumsum(rev(lv)))
  expect_equal(cslv, unname(cslv2[-1]))
})



# Create Data useful for many tests
set.seed(1)
samples <- rnbinom(n = 30, mu = 4, size = 1/0.7)
mu <- rnorm(n = 30, mean = 4)
X <- matrix(rnorm(n = 30 * 4), nrow = 30, ncol = 4)


test_that("Convential and Bandara approach produce the same estimates", {
  pbl_w_cr <- gampoi_overdispersion_mle(y = samples, mean = mu,
                                        model_matrix = X, do_cox_reid_adjustment = TRUE)$estimate
  ban_w_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = TRUE)$estimate
  con_w_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = TRUE)$estimate
  expect_equal(ban_w_cr, con_w_cr, tolerance = 1e-4)
  expect_equal(ban_w_cr, pbl_w_cr, tolerance = 1e-4)

  pbl_wo_cr <- gampoi_overdispersion_mle(y = samples, mean = mu,
                                        model_matrix = X, do_cox_reid_adjustment = FALSE)$estimate
  ban_wo_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = FALSE)$estimate
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = FALSE)$estimate
  expect_equal(ban_wo_cr, con_wo_cr, tolerance = 1e-4)
  expect_equal(ban_wo_cr, pbl_wo_cr, tolerance = 1e-4)
})


test_that("Derivative of Bandara score function works", {
  rg <- seq(1, 30, length.out = 1001)
  score_values <- sapply(rg, function(r){
    score_function_bandara_fast(samples, makeCumSumLookupVector(samples), mu = mu, r = r,
                                model_matrix = X, do_cr_adj = TRUE)
  })
  emp_deriv<- diff(score_values) / diff(rg)[1]


  analyt_deriv <- sapply(rg, function(r){
    score_deriv_function_bandara_fast(samples, makeCumSumLookupVector(samples), mu = mu, r = r,
                                      model_matrix = X, do_cr_adj = TRUE)
  })
  respace_analyt <- zoo::rollmean(analyt_deriv, 2)
  fit <- lm(respace_analyt ~ emp_deriv)
  expect_equal(unname(coef(fit)["(Intercept)"]), 0, tolerance = 1e-4)
  expect_equal(unname(coef(fit)[2]), 1, tolerance = 1e-3)
})



test_that("C++ implementation of loglikelihood and score match", {
  log_theta_g <- seq(-3, 3, length.out = 1001)

  ll_values <- sapply(log_theta_g, function(log_theta){
    conventional_loglikelihood_fast(samples, mu = mu, log_theta = log_theta, model_matrix = X, do_cr_adj = TRUE)
  })
  emp_deriv<- diff(ll_values) / diff(log_theta_g)[1]
  analyt_deriv <- sapply(log_theta_g, function(log_theta){
    conventional_score_function_fast(samples, mu = mu, log_theta = log_theta, model_matrix = X, do_cr_adj = TRUE)
  })
  respace_analyt <- zoo::rollmean(analyt_deriv, 2)
  expect_equal(emp_deriv, respace_analyt, tolerance = 1e-5)
})


test_that("C++ implementation of score and score_deriv match", {
  log_theta_g <- seq(-3, 3, length.out = 1001)

  score_values <- sapply(log_theta_g, function(log_theta){
    conventional_score_function_fast(samples, mu = mu, log_theta = log_theta, model_matrix = X, do_cr_adj = TRUE)
  })
  emp_deriv<- diff(score_values) / diff(log_theta_g)[1]
  analyt_deriv <- sapply(log_theta_g, function(log_theta){
    conventional_deriv_score_function_fast(samples, mu = mu, log_theta = log_theta, model_matrix = X, do_cr_adj = TRUE)
  })
  respace_analyt <- zoo::rollmean(analyt_deriv, 2)
  expect_equal(emp_deriv, respace_analyt, tolerance = 1e-5)
})




test_that("Estimation methods can handle under-dispersion", {

  # I know that for seed 1, this produces an under-dispersed sample
  set.seed(1)
  samples <- rpois(n = 100, lambda = 5)
  expect_gt(mean(samples), var(samples))  # Proof of underdispersion
  mu <- rep(mean(samples), length(samples))
  ban_wo_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)$estimate
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)$estimate

  expect_equal(ban_wo_cr, 0)
  expect_equal(con_wo_cr, 0)

  # However, if mu is large then again a theta can be estimated
  mu <- rep(10, length(samples))
  ban_wo_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                          do_cox_reid_adjustment = FALSE)$estimate
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                               do_cox_reid_adjustment = FALSE)$estimate
  expect_true(ban_wo_cr != 0 && con_wo_cr != 0)
  expect_equal(ban_wo_cr, con_wo_cr)
})


test_that("Estimation methods can handle Infinite dispersion", {
  # For some reason this model matrix makes the dispersion estimate
  # go to +Inf. Fixed by adding a cr_correction_factor = 0.99 to
  # the calculation of the Cox-Reid adjustment.
  # The problem was that the lgamma(1/theta) and CR-adjustment
  # canceled each other exactly
  model_matrix <- cbind(1, rnorm(n=5))
  mean_vector <- c(0.2, 0.6, 0.8, 0.2, 0.1)
  y <- c(0, 0, 3, 0, 0)

  fit <- gampoi_overdispersion_mle(y, mean_vector, model_matrix = model_matrix)
  expect_lt(fit$estimate, 1e5)
})


test_that("Estimation methods can handle mu = 0", {

  mu <- c(head(mu, length(samples)-1), 0)

  b_res <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                             do_cox_reid_adjustment = FALSE)
  c_res <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                  do_cox_reid_adjustment = FALSE)
  expect_equal(b_res$estimate, c_res$estimate)

})


test_that("Identical y values work", {
  # Underdispersed data -> theta = 0
  samples <- rep(6, times = 10)
  mu <- rep(mean(samples), length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                             do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                  do_cox_reid_adjustment = FALSE)
  expect_equal(ban$estimate, 0)
  expect_equal(con$estimate, 0)

  # If mu is small enough, it works again
  samples <- rep(6, times = 10)
  mu <- rep(1, length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                    do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)
  expect_equal(ban$estimate, con$estimate, tolerance = 1e-4)


  # For all y = 0 -> cannot really make inference, assume theta = 0
  samples <- rep(0, times = 10)
  mu <- rep(1, length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                    do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)
  expect_equal(ban$estimate, 0)
  expect_equal(con$estimate, 0)

})

test_that("one value is enough to get an answer", {
  expect_equal(gampoi_overdispersion_mle(y = 3)$estimate, 0)
  expect_equal(bandara_overdispersion_mle(y = 3, mean_vector = 3)$estimate, 0)
  expect_equal(conventional_overdispersion_mle(y = 3, mean_vector = 3)$estimate, 0)
  expect_false(bandara_overdispersion_mle(y = 3, mean_vector = 1.3)$estimate == 0)
  expect_equal(bandara_overdispersion_mle(y = 3, mean_vector = 1.3)$estimate,
               conventional_overdispersion_mle(y = 3, mean_vector = 1.3)$estimate, tolerance = 1e-6)
})


test_that("subsampling works and does not affect performance too badly", {
  y <- rnbinom(n = 1e4, mu = 5, size = 1 / 0.7)
  r1 <- gampoi_overdispersion_mle(y = y, n_subsamples = 1e4)
  r2 <- gampoi_overdispersion_mle(y = y, n_subsamples = 1000)
  expect_lt(abs(r1$estimate - r2$estimate), 0.1)

})


test_that("weird subsampling values are handled correctly", {
  y <- rnbinom(n = 20, mu = 5, size = 1 / 0.7)
  r1 <- gampoi_overdispersion_mle(y = y, n_subsamples = 1)
  r0 <- gampoi_overdispersion_mle(y = y, n_subsamples = 0)
  expect_equal(r0$estimate, 0)
  expect_error(
    gampoi_overdispersion_mle(y = y, n_subsamples = -3)
  )
  expect_error(
    gampoi_overdispersion_mle(y = y, n_subsamples = c(3, 10))
  )
  rfull <- gampoi_overdispersion_mle(y = y)
  r3000 <- gampoi_overdispersion_mle(y = y, n_subsamples = 3000)
  rInf <- gampoi_overdispersion_mle(y = y, n_subsamples = Inf)
  expect_equal(rfull, r3000)
  expect_equal(rfull, rInf)
})



test_that("DelayedArrays are handled efficiently", {
  n_genes <- 100
  n_samples <- 40
  model_matrix <- matrix(1, nrow = n_samples, ncol = 1)
  mat <- matrix(seq_len(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
  mat_hdf5 <-  as(mat, "HDF5Matrix")
  offset_matrix <- combine_size_factors_and_offset(TRUE, 1, mat_hdf5)$offset_matrix
  dispersions <- estimate_dispersions_roughly(mat_hdf5, model_matrix, offset_matrix)

  beta_vec_init <- estimate_betas_roughly_one_group(mat_hdf5, offset_matrix)
  Betas <- estimate_betas_one_group(mat_hdf5, offset_matrix, dispersions, beta_vec_init)$Beta
  mean_matrix <- calculate_mu(Betas, model_matrix, offset_matrix)
  mean_matrix_ram <- as.matrix(mean_matrix)


  disp_est_r_ram <- vapply(seq_len(n_genes), function(gene_idx){
    gampoi_overdispersion_mle(y = mat_hdf5[gene_idx, ], mean = mean_matrix[gene_idx, ],
                              model_matrix = model_matrix, do_cox_reid_adjustment = TRUE,
                              n_subsamples = n_samples)$estimate
  }, FUN.VALUE = 0.0)
  disp_est_r_hdf5 <- vapply(seq_len(n_genes), function(gene_idx){
    gampoi_overdispersion_mle(y = mat_hdf5[gene_idx, ], mean = mean_matrix[gene_idx, ],
                              model_matrix = model_matrix, do_cox_reid_adjustment = TRUE,
                              n_subsamples = n_samples)$estimate
  }, FUN.VALUE = 0.0)

  beachmat_ram <- estimate_overdispersions_fast(mat, mean_matrix_ram, model_matrix,
                                  do_cox_reid_adjustment = TRUE, n_subsamples = n_samples)
  beachmat_hdf5 <- estimate_overdispersions_fast(mat_hdf5, mean_matrix, model_matrix,
                                  do_cox_reid_adjustment = TRUE, n_subsamples = n_samples)

  expect_equal(disp_est_r_ram, disp_est_r_hdf5)
  expect_equal(disp_est_r_ram, beachmat_ram)
  expect_equal(disp_est_r_ram, beachmat_hdf5)
})







