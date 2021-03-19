test_that("residual calculation works", {
  set.seed(1)
  X <- cbind(1, matrix(rnorm(4 * 2) , nrow = 4, ncol = 2))
  Y <- matrix(rnbinom(n = 2 * 4, mu = 30, size = 0.7), nrow = 2, ncol = 4)
  fit <- glm_gp(Y, X, size_factors = FALSE, overdispersion = 1/0.7)
  r_fit1 <- glm(Y[1,] ~ X - 1, family = MASS::negative.binomial(theta = 0.7))
  r_fit2 <- glm(Y[2,] ~ X - 1, family = MASS::negative.binomial(theta = 0.7))
  expect_equal(unname(fit$Beta[1,]), unname(coef(r_fit1)), tolerance = 1e-4)
  expect_equal(unname(fit$Beta[2,]), unname(coef(r_fit2)), tolerance = 1e-4)

  expect_equal(c(t(residuals(fit, "response"))),
              unname(c(residuals.glm(r_fit1, "response"), residuals.glm(r_fit2, "response"))),
              tolerance = 1e-4)
  expect_equal(c(t(residuals(fit, "working"))),
               unname(c(residuals.glm(r_fit1, "working"), residuals.glm(r_fit2, "working"))),
               tolerance = 1e-5)
  expect_equal(c(t(residuals(fit, "pearson"))),
               unname(c(residuals.glm(r_fit1, "pearson"), residuals.glm(r_fit2, "pearson"))),
               tolerance = 1e-5)
  expect_equal(c(t(residuals(fit, "deviance"))),
               unname(c(residuals.glm(r_fit1, "deviance"), residuals.glm(r_fit2, "deviance"))),
               tolerance = 1e-5)

  # Randomized Quantiles are by definition not equal
  r_qs <- c(statmod::qresiduals(r_fit1),  statmod::qresiduals(r_fit2))
  res <- c(t(residuals(fit, "randomized_quantile")))
  expect_gt(cor(r_qs, res), 0.99)
})


test_that("residuals are named", {
  set.seed(1)
  X <- cbind(1, matrix(rnorm(4 * 2) , nrow = 4, ncol = 2))
  Y <- matrix(rnbinom(n = 2 * 4, mu = 30, size = 0.7), nrow = 2, ncol = 4)
  rownames(Y) <- paste0("Gene_", seq_len(nrow(Y)))
  colnames(Y) <- paste0("Cell_", seq_len(ncol(Y)))
  fit <- glm_gp(Y, X, size_factors = FALSE, overdispersion = 1/0.7)
  expect_equal(dimnames(residuals(fit)), dimnames(Y))
})



test_that("residual calculation works with Delayed Matrix", {
  set.seed(1)
  X <- cbind(1, matrix(rnorm(4 * 2) , nrow = 4, ncol = 2))
  Y <- matrix(rnbinom(n = 2 * 4, mu = 30, size = 0.7), nrow = 2, ncol = 4)
  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)
  fit <- glm_gp(Y_hdf5, X, size_factors = FALSE, overdispersion = 1/0.7)
  r_fit1 <- glm(Y[1,] ~ X - 1, family = MASS::negative.binomial(theta = 0.7))
  r_fit2 <- glm(Y[2,] ~ X - 1, family = MASS::negative.binomial(theta = 0.7))
  expect_equal(unname(fit$Beta[1,]), unname(coef(r_fit1)), tolerance = 1e-4)
  expect_equal(unname(fit$Beta[2,]), unname(coef(r_fit2)), tolerance = 1e-4)

  expect_s4_class(residuals(fit, "response"), "DelayedMatrix")
  expect_s4_class(residuals(fit, "working"), "DelayedMatrix")
  expect_s4_class(residuals(fit, "pearson"), "DelayedMatrix")
  expect_s4_class(residuals(fit, "deviance"), "DelayedMatrix")

  expect_true(DelayedArray::isPristine(residuals(fit, "response")))
  expect_true(DelayedArray::isPristine(residuals(fit, "working")))
  expect_true(DelayedArray::isPristine(residuals(fit, "pearson")))
  expect_true(DelayedArray::isPristine(residuals(fit, "deviance")))


  expect_equal(c(t(residuals(fit, "response"))),
               unname(c(residuals.glm(r_fit1, "response"), residuals.glm(r_fit2, "response"))),
               tolerance = 1e-4)
  expect_equal(c(t(residuals(fit, "working"))),
               unname(c(residuals.glm(r_fit1, "working"), residuals.glm(r_fit2, "working"))),
               tolerance = 1e-5)
  expect_equal(c(t(residuals(fit, "pearson"))),
               unname(c(residuals.glm(r_fit1, "pearson"), residuals.glm(r_fit2, "pearson"))),
               tolerance = 1e-5)
  expect_equal(c(t(residuals(fit, "deviance"))),
               unname(c(residuals.glm(r_fit1, "deviance"), residuals.glm(r_fit2, "deviance"))),
               tolerance = 1e-5)

  # Randomized Quantiles are by definition not equal
  r_qs <- c(statmod::qresiduals(r_fit1),  statmod::qresiduals(r_fit2))
  res <- c(t(residuals(fit, "randomized_quantile")))
  expect_gt(cor(r_qs, res), 0.99)

})




test_that("qres.gampoi can handle extreme value", {
  Y <- matrix(27)
  Mu <- matrix(2)
  overdispersion <- 0
  res <- qres.gampoi(Y, Mu, overdispersion)
  expect_false(is.infinite(res))

  Y <- matrix(2700)
  res <- qres.gampoi(Y, Mu, overdispersion)
  expect_false(is.infinite(res))


  Y <- matrix(2)
  Mu <- matrix(270)
  res <- qres.gampoi(Y, Mu, overdispersion)
  expect_false(is.infinite(res))


  Y <- matrix(c(2, 2), ncol = 1)
  Mu <- matrix(c(270, 270), ncol = 1)
  overdispersion <- c(0, 500)
  res <- qres.gampoi(Y, Mu, overdispersion)
  expect_false(is.infinite(res[1,1]), is.infinite(res[2,1]))
  expect_false(res[1,1] == res[2,1])
})



test_that("qres.gampoi can handle other weird values", {
  # This specific combination of parameters caused NA's
  Y <- matrix(27)
  Mu <- matrix(0.435023)
  overdispersion <- 0

  a <- ppois(Y - 1, lambda = Mu)
  b <- ppois(Y, lambda = Mu)
  # This really shouldn't happen
  # Nonetheless, it does for this combination of parameters
  # at least on Linux and MacOS
  if(! is_windows()){
    expect_gt(a, b)
  }

  # However, now qres.gampoi handles this edge case
  res <- qres.gampoi(Y, Mu, overdispersion)
  expect_false(is.nan(res))
})



test_that("compute_gp_deviance can handle weird values", {
  # This specific combination of parameters caused negative deviance
  y <- 1
  mu <- 0.99999999999994
  theta <- 1e-7

  # However, now qres.gampoi handles this edge case
  res <- compute_gp_deviance(y, mu, theta)
  expect_gte(res, 0)
})


