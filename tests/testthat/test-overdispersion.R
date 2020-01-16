


test_that("Overdispersion estimator works", {
  samples <- distraltparam::raltnbinom(n = 100, mean = 4, dispersion = 0.7)
  mu <- rep(mean(samples), length(samples))
  X <- matrix(1, nrow = length(samples), ncol = 1)
  res1 <- gampoi_overdispersion_mle(y = samples, mean_vector = mu)$root
  res2 <- optimize(function(r){
    -  (sum(distraltparam::daltnbinom(x = samples, mean=mu, dispersion = 1/r, log=TRUE)) +
          0.5 * log(det(t(X) %*% diag(1/(1/mu + r)) %*% X)))
  }, lower=1e-3, upper=1000)$minimum^(-1)
  expect_equal(res1, res2, tolerance = 1e-3)
})
