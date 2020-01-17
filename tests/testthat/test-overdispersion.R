


test_that("makeCumSumLookupVector works", {
  vector <- rpois(n = 100, lambda = 10)
  cslv <- makeCumSumLookupVector(vector)
  lv <- table(factor(vector, levels = seq_len(max(vector)+1)-1, ordered = TRUE))
  cslv2 <- rev(cumsum(rev(lv)))
  expect_equal(cslv, unname(cslv2[-1]))
})



# Create Data useful for many tests

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




