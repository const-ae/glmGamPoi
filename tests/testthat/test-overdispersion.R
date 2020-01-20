


test_that("makeCumSumLookupVector works", {
  vector <- rpois(n = 100, lambda = 10)
  cslv <- makeCumSumLookupVector(vector)
  lv <- table(factor(vector, levels = seq_len(max(vector)+1)-1, ordered = TRUE))
  cslv2 <- rev(cumsum(rev(lv)))
  expect_equal(cslv, unname(cslv2[-1]))
})



# Create Data useful for many tests
set.seed(1)
samples <- distraltparam::raltnbinom(n = 30, mean = 4, dispersion = 0.7)
mu <- rnorm(n = 30, mean = 4)
X <- matrix(rnorm(n = 30 * 4), nrow = 30, ncol = 4)


test_that("Convential and Bandara approach produce the same estimates", {
  ban_w_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = TRUE)$root
  con_w_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = TRUE)$root
  expect_equal(ban_w_cr, con_w_cr, tolerance = 1e-4)

  ban_wo_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = FALSE)$root
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = FALSE)$root
  expect_equal(ban_wo_cr, con_wo_cr, tolerance = 1e-4)
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
                                         do_cox_reid_adjustment = FALSE)$root
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)$root

  expect_equal(ban_wo_cr, 0)
  expect_equal(con_wo_cr, 0)

  # However, if mu is large then again a theta can be estimated
  mu <- rep(10, length(samples))
  ban_wo_cr <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                          do_cox_reid_adjustment = FALSE)$root
  con_wo_cr <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                               do_cox_reid_adjustment = FALSE)$root
  expect_true(ban_wo_cr != 0 && con_wo_cr != 0)
  expect_equal(ban_wo_cr, con_wo_cr)
})



test_that("Estimation methods can handle mu = 0", {

  mu <- c(head(mu, length(samples)-1), 0)

  expect_error({
    bandara_overdispersion_mle(y = samples, mean_vector = mu,
                               do_cox_reid_adjustment = FALSE)
  })
  expect_error({
    conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                    do_cox_reid_adjustment = FALSE)
  })

})


test_that("Identical y values work", {
  # Underdispersed data -> theta = 0
  samples <- rep(6, times = 10)
  mu <- rep(mean(samples), length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                             do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                  do_cox_reid_adjustment = FALSE)
  expect_equal(ban$root, 0)
  expect_equal(con$root, 0)

  # If mu is small enough, it works again
  samples <- rep(6, times = 10)
  mu <- rep(1, length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                    do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)
  expect_equal(ban$root, con$root, tolerance = 1e-4)


  # For all y = 0 -> cannot really make inference, assume theta = 0
  samples <- rep(0, times = 10)
  mu <- rep(1, length(samples))
  ban <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                    do_cox_reid_adjustment = FALSE)
  con <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                         do_cox_reid_adjustment = FALSE)
  expect_equal(ban$root, 0)
  expect_equal(con$root, 0)

})


