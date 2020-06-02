


test_that("make_table works", {
  vector <- rpois(n = 100, lambda = 10)
  tab <- make_table(vector)
  tab2 <- table(vector)
  expect_equal(as.numeric(tab2), tab[[2]][order(tab[[1]])])
  expect_equal(as.numeric(names(tab2)), sort(tab[[1]]))
})


test_that("digamma approximation works", {
  # In the derivative of the loglikelihood wrt to log theta
  # the term (sum(digamma(y + x))  - length(y) * digamma(x)) * x appears
  # This is always smaller than sum(y) and for large values of x it is
  # approximately equal to sum(y).
  x <- 1e6
  y <- c(3, 1, 6)
  expect_equal((sum(digamma(y + x))  - length(y) * digamma(x)) * x, sum(y), tolerance = 1e-4)
  x <- 1e-5
  expect_lt((sum(digamma(y + x))  - length(y) * digamma(x)) * x, sum(y))

  ## Due to numerical imprecision at very large numbers it can look as if
  ## the left term would become larger, but that is wrong.
  # x <- 1e15
  # expect_lt((sum(digamma(y + x))  - length(y) * digamma(x)) * x, sum(y))
  ## Check out the plot
  # xg <- seq(-35, 35, l = 1001)
  # values <- sapply(exp(xg), function(x) (sum(digamma(y + x))  - length(y) * digamma(x)) * x)
  # plot(xg, values, col = (values < sum(y)) + 1)
  # abline(h = sum(y))
})


# Create Data useful for many tests
set.seed(1)
samples <- rnbinom(n = 30, mu = 4, size = 1/0.7)
mu <- rnorm(n = 30, mean = 4)
X <- matrix(rnorm(n = 30 * 4), nrow = 30, ncol = 4)



test_that("Score function can handle extreme inputs properly", {
  samples <- rpois(n = 3, lambda = 3)
  mu <- rnorm(n = length(samples), mean = 4)
  X <- matrix(rnorm(n = length(samples) * 4), nrow = length(samples), ncol = 4)
  tab <- make_table(samples)
  expect_equal(conventional_score_function_fast(samples, mu, -35, X, do_cr_adj = TRUE),
               conventional_score_function_fast(samples, mu, -35, X, do_cr_adj = TRUE, tab[[1]], tab[[2]]))

  # xg <- seq(-30, 10, l = 51)
  # values1 <- sapply(xg, function(lt) conventional_score_function_fast2(samples, mu, lt, X[,1,drop=FALSE], do_cr_adj = FALSE, tab[[1]], tab[[2]]))
  # values2 <- sapply(xg, function(lt) conventional_score_function_fast(samples, mu, lt, X[,1,drop=FALSE], do_cr_adj = FALSE, tab[[1]], tab[[2]]))
  # plot(xg, values1)
  # lines(xg, values2, col = "red")
  # plot(xg, values2, col = (values2 < 0) + 1)

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
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)$estimate

  expect_equal(con_wo_cr, 0)

  # However, if mu is large then again a theta can be estimated
  mu <- rep(10, length(samples))
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                               do_cox_reid_adjustment = FALSE)$estimate
  expect_true(con_wo_cr != 0)
})


test_that("gampoi_overdispersion_mle can handle weird input 1", {
  y <- c(5, 2, 1, 1, 3, 1, 1)
  X <- cbind(1, c(0, 0, 0, 1, 1, 1, 1))
  mu <- c(2.6, 2.6, 2.6, 1.5, 1.5, 1.5, 1.5)
  # This used to fail  because no starting position is found
  # because mean was exactly equal to var
  est_2 <- conventional_overdispersion_mle(y, mean_vector = mu,
                                  model_matrix = X,
                                  do_cox_reid_adjustment = TRUE,
                                  verbose = FALSE)
  expect_true(TRUE)
})


test_that("gampoi_overdispersion_mle can handle weird input 2", {
  y <- c(5, 10, 3, 4, 5, 4, 5)
  X <- cbind(1, c(0, 0, 0, 1, 1, 1, 1))
  mu <- c(6, 6, 6, 4.5, 4.5, 4.5, 4.5)
  # This used to fail because mean was exactly equal to var
  est_2 <- conventional_overdispersion_mle(y, mean_vector = mu,
                                           model_matrix = X,
                                           do_cox_reid_adjustment = TRUE,
                                           verbose = FALSE)
  expect_true(TRUE)
})

test_that("gampoi_overdispersion_mle can handle weird input 3", {
  y <- c(rep(0, times = 399), 10)
  X <- cbind(1, sample(c(0,1), 400, replace = TRUE), rnorm(400))
  mu <- rep(0.7, 400)
  est <- conventional_overdispersion_mle(y, mean_vector = mu,
                                           model_matrix = X,
                                           do_cox_reid_adjustment = TRUE,
                                           verbose = FALSE)

  expect_gt(est$estimate, 1e8)
})

test_that("gampoi_overdispersion_mle can handle weird input 4", {
  y <- c(rep(0, times = 5), 2)
  mu <- c(rep(1e-30, 5), 1.999999999999976)
  X <- cbind(1, rep(c(0,1), each = 3), rnorm(6))
  log_theta <- -3
  res <- conventional_loglikelihood_fast(y, mu = mu, log_theta = log_theta,
                                  model_matrix = X, do_cr_adj = TRUE)

  expect_true(is.finite(res))

  score <- conventional_score_function_fast(y, mu = mu, log_theta = log_theta,
                                  model_matrix = X, do_cr_adj = TRUE)
  expect_true(is.finite(score))

  deriv <- conventional_deriv_score_function_fast(y, mu = mu, log_theta = log_theta,
                                            model_matrix = X, do_cr_adj = TRUE)
  expect_true(is.finite(deriv))
  # w <- 1/(1/mu + exp(log_theta) + 1e-6)
  # b <- t(X) %*% diag(w) %*% X
  # det(b)
  # xg <- seq(-25, 25, l = 1001)
  # values <- sapply(xg, function(lt){
  #   conventional_score_function_fast(y, mu = mu, log_theta = lt,
  #                                   model_matrix = X, do_cr_adj = TRUE)
  # })
  # plot(xg, values)
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

  c_res <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                  do_cox_reid_adjustment = FALSE)
  expect_true(TRUE)
})


test_that("Identical y values work", {
  # Underdispersed data -> theta = 0
  samples <- rep(6, times = 10)
  mu <- rep(mean(samples), length(samples))
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                  do_cox_reid_adjustment = FALSE)
  expect_equal(con$estimate, 0)

  # If mu is small enough, it works again
  samples <- rep(6, times = 10)
  mu <- rep(1, length(samples))
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)


  # For all y = 0 -> cannot really make inference, assume theta = 0
  samples <- rep(0, times = 10)
  mu <- rep(1, length(samples))
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)
  expect_equal(con$estimate, 0)

})

test_that("one value is enough to get an answer", {
  expect_equal(gampoi_overdispersion_mle(y = 3, mean = 3.0001)$estimate, 0)
  expect_equal(conventional_overdispersion_mle(y = 3, mean_vector = 3.0001)$estimate, 0)
})


test_that("subsampling works and does not affect performance too badly", {
  y <- rnbinom(n = 1e4, mu = 5, size = 1 / 0.7)
  r1 <- gampoi_overdispersion_mle(y = y, subsample = 1e4)
  r2 <- gampoi_overdispersion_mle(y = y, subsample = 1000)
  expect_lt(abs(r1$estimate - r2$estimate), 0.1)

})


test_that("weird subsampling values are handled correctly", {
  y <- rnbinom(n = 20, mu = 5, size = 1 / 0.7)
  r1 <- gampoi_overdispersion_mle(y = y, subsample = 1)
  r0 <- gampoi_overdispersion_mle(y = y, subsample = 0)
  expect_equal(r0$estimate, 0)
  expect_error(
    gampoi_overdispersion_mle(y = y, subsample = -3)
  )
  expect_error(
    gampoi_overdispersion_mle(y = y, subsample = c(3, 10))
  )
  rfull <- gampoi_overdispersion_mle(y = y)
  r3000 <- gampoi_overdispersion_mle(y = y, subsample = 3000)
  rInf <- gampoi_overdispersion_mle(y = y, subsample = Inf)
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
                              model_matrix = model_matrix, do_cox_reid_adjustment = TRUE)$estimate
  }, FUN.VALUE = 0.0)
  disp_est_r_hdf5 <- vapply(seq_len(n_genes), function(gene_idx){
    gampoi_overdispersion_mle(y = mat_hdf5[gene_idx, ], mean = mean_matrix[gene_idx, ],
                              model_matrix = model_matrix, do_cox_reid_adjustment = TRUE)$estimate
  }, FUN.VALUE = 0.0)

  beachmat_ram <- estimate_overdispersions_fast(mat, mean_matrix_ram, model_matrix,
                                  do_cox_reid_adjustment = TRUE, n_subsamples = n_samples)
  beachmat_hdf5 <- estimate_overdispersions_fast(mat_hdf5, mean_matrix, model_matrix,
                                  do_cox_reid_adjustment = TRUE, n_subsamples = n_samples)

  expect_equal(disp_est_r_ram, disp_est_r_hdf5)
  expect_equal(disp_est_r_ram, beachmat_ram)
  expect_equal(disp_est_r_ram, beachmat_hdf5)
})







