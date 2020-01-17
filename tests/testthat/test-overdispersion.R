


test_that("makeCumSumLookupVector works", {
  vector <- rpois(n = 100, lambda = 10)
  cslv <- makeCumSumLookupVector(vector)
  lv <- table(factor(vector, levels = seq_len(max(vector)+1)-1, ordered = TRUE))
  cslv2 <- rev(cumsum(rev(lv)))
  expect_equal(cslv, unname(cslv2[-1]))
})


test_that("Convential and Bandara approach produce the same estimates", {
  samples <- distraltparam::raltnbinom(n = 30, mean = 4, dispersion = 0.7)
  mu <- rep(mean(samples), length(samples))
  X <- matrix(1, nrow = length(samples), ncol = 1)

  res1 <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = TRUE)$root
  res2 <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = TRUE)$root
  expect_equal(res1, res2, tolerance = 1e-4)

  res1 <- bandara_overdispersion_mle(y = samples, mean_vector = mu,
                                     model_matrix = X, do_cox_reid_adjustment = FALSE)$root
  res2 <- conventional_overdispersion_mle(y = samples, mean_vector = mu,
                                          model_matrix = X, do_cox_reid_adjustment = FALSE)$root
  expect_equal(res1, res2, tolerance = 1e-4)
})



test_that("Derivative of Bandara score function works", {
  samples <- distraltparam::raltnbinom(n = 30, mean = 4, dispersion = 0.7)
  mu <- rep(mean(samples), length(samples))
  X <- matrix(1, nrow = length(samples), ncol = 1)

  rg <- seq(3, 30, l=1001)
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



