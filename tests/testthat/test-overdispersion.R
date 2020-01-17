


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

