

test_that("loc_median_fit works", {
  x <- sample(1:10)
  y <- rnorm(n = 10)
  r1 <- loc_median_fit(x, y)

  r2 <- loc_median_fit(x, y, fraction = 1, weighted = FALSE)
  expect_equal(r2, rep(median(y), 10))
  r25 <- loc_median_fit(x[-1], y[-1], fraction = 1, weighted = FALSE)
  expect_equal(r25, rep(median(y[-1]), 9))

  r4 <- loc_median_fit(x[1], y[1])
  expect_equal(r4, y[1])

  r5 <- loc_median_fit(numeric(0), numeric(0))
  expect_equal(r5, numeric(0))
})



test_that("shrink_ql_dispersion is robust", {

  expect_silent({
    mu <- rnorm(n = 100, mean = 5, sd = 0.1)
    disp_est <- rnorm(n = 100,  mean = 2, sd = 0.2)
    my_res <- shrink_ql_dispersion(disp_est, mu, df = 2)


    mu <- rnorm(n = 100, mean = 5, sd = 0.1)
    disp_est <- rep(2.1, 100)
    shrink_ql_dispersion(disp_est, mu, df = 2)


    mu <- c(3, 2)
    disp_est <- c(1,1)
    shrink_ql_dispersion(disp_est, mu, df = 2)

    mu <- c(3, 0)
    disp_est <- c(1,1)
    shrink_ql_dispersion(disp_est, mu, df = 2)


    mu <- 2
    disp_est <- 2
    shrink_ql_dispersion(disp_est, mu, df = 2)

    mu <- 0
    disp_est <- 1
    shrink_ql_dispersion(disp_est, mu, df = 2)
  })
})


test_that("variance prior estimation works", {
  n <- 1000
  true_df0 <- 5
  true_variance0 <- 1.3
  # Based on the first equation of section 3 in the 2004 limma paper by Smyth
  true_variance <- 1/(rchisq(n = n, df = true_df0) / (true_df0 * true_variance0))
  observations <- t(vapply(seq_len(n), function(idx){
    rnorm(n = 5, mean = 0, sd = sqrt(true_variance[idx]))
  }, FUN.VALUE = rep(0.0, 5)))

  obs_var <- DelayedMatrixStats::rowVars(observations)
  df <- ncol(observations) - 1
  res_limma <- limma:::squeezeVar(obs_var, df = df)
  res_gp <- variance_prior(obs_var, df)

  loglikelihood_gp <- sum(df(obs_var/res_gp$variance0, df1=df, df2=res_gp$df0, log=TRUE) - log(res_gp$variance0))
  loglikelihood_limma <- sum(df(obs_var/res_limma$var.prior, df1=df, df2=res_limma$df.prior, log=TRUE) - log(res_limma$var.prior))
  expect_gt(loglikelihood_gp, loglikelihood_limma)

})

test_that("variance prior estimation works with covariates", {
  n <- 1000
  true_df0 <- 50
  covariate <- seq_len(n) / n
  true_variances0 <- 20 * covariate + 0.3
  # Based on the first equation of section 3 in the 2004 limma paper by Smyth
  true_variance <- 1/(rchisq(n = n, df = true_df0) / (true_df0 * true_variances0))
  observations <- t(vapply(seq_len(n), function(idx){
    rnorm(n = 5, mean = 0, sd = sqrt(true_variance[idx]))
  }, FUN.VALUE = rep(0.0, 5)))

  obs_var <- DelayedMatrixStats::rowVars(observations)
  df <- ncol(observations) - 1
  res_limma <- limma::squeezeVar(obs_var, df = df, covariate = covariate)
  res_gp <- variance_prior(obs_var, df, covariate = covariate)
  res_gp2 <- variance_prior(obs_var, df, covariate = NULL)


  loglikelihood_gp <- sum(df(obs_var/res_gp$variance0, df1=df, df2=res_gp$df0, log=TRUE) - log(res_gp$variance0))
  loglikelihood_limma <- sum(df(obs_var/res_limma$var.prior, df1=df, df2=res_limma$df.prior, log=TRUE) - log(res_limma$var.prior))
  expect_gt(loglikelihood_gp, loglikelihood_limma)

})



test_that("gampoi_test_qlr works", {
  set.seed(1)
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  fit <- glm_gp(Y, design = design)
  res <- gampoi_test_qlr(Y, fit, reduced_design = ~ group + cont1, col_data = annot)
  res2 <- gampoi_test_qlr(Y, fit, contrast = cont2)
  # Should be equal except for log-fold change column
  expect_equal(res[,-7], res2[,-7], tolerance = 1e-6)
  expect_equal(res$lfc, rep(NA, 30), tolerance = 1e-6)
  expect_equal(res2$lfc, fit$Beta[,"cont2"] / log(2), tolerance = 1e-6)

  design_wo_intercept <- model.matrix(~ group + cont1 + cont2 - 1, data = annot)
  fit_wo_intercept <- glm_gp(Y, design = design_wo_intercept)
  res3 <- gampoi_test_qlr(Y, fit_wo_intercept, reduced_design = ~ cont1 + cont2 + 1, col_data = annot)
  res4 <- gampoi_test_qlr(Y, fit_wo_intercept, contrast = groupA - groupB)
  expect_equal(res3[,-7], res4[,-7], tolerance = 1e-6)

  res5 <- gampoi_test_qlr(Y, fit, contrast = "-groupB", col_data = annot)
  expect_equal(res4, res5, tolerance = 0.1)

})


test_that("gampoi_test_qlr works with contrast vector input", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  fit <- glm_gp(Y, design = design)
  res1 <- gampoi_test_qlr(Y, fit, contrast = cont2 - groupB)
  res2 <- gampoi_test_qlr(Y, fit, contrast = c(0, -1, 0, 1))
  expect_equal(res1, res2)

})





