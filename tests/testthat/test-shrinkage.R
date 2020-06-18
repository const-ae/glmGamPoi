

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



test_that("overdispersion_shrinkage is robust", {

  expect_silent({
    mu <- rnorm(n = 100, mean = 5, sd = 0.1)
    disp_est <- rnorm(n = 100,  mean = 2, sd = 0.2)
    my_res <- overdispersion_shrinkage(disp_est, mu, df = 2)


    mu <- rnorm(n = 100, mean = 5, sd = 0.1)
    disp_est <- rep(2.1, 100)
    overdispersion_shrinkage(disp_est, mu, df = 2)


    mu <- c(3, 2)
    disp_est <- c(1,1)
    overdispersion_shrinkage(disp_est, mu, df = 2)

    mu <- c(3, 0)
    disp_est <- c(1,1)
    overdispersion_shrinkage(disp_est, mu, df = 2)


    mu <- 2
    disp_est <- 2
    overdispersion_shrinkage(disp_est, mu, df = 2)

    mu <- 0
    disp_est <- 1
    overdispersion_shrinkage(disp_est, mu, df = 2)
  })
})


test_that("spline fit is not degenerate", {

  log_covariate <- sample(c(rep(-25, 350), rnorm(n = 650, mean = 0, sd = 2.6)))
  ra <- range(log_covariate)
  knots <- ra[1] + c(1/3, 2/3) * (ra[2] - ra[1])
  design <- splines::ns(log_covariate, df = 4, knots = knots, intercept = TRUE)
  expect_equal(qr(design)$rank, 4)

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
  res_gp <- variance_prior(obs_var, df, covariate = covariate, abundance_trend = TRUE)
  res_gp2 <- variance_prior(obs_var, df, covariate = NULL)


  loglikelihood_gp <- sum(df(obs_var/res_gp$variance0, df1=df, df2=res_gp$df0, log=TRUE) - log(res_gp$variance0))
  loglikelihood_limma <- sum(df(obs_var/res_limma$var.prior, df1=df, df2=res_limma$df.prior, log=TRUE) - log(res_limma$var.prior))
  expect_gt(loglikelihood_gp, loglikelihood_limma)

})

test_that("variance prior estimation works with zeros", {

  n <- 1000 + 10
  true_df0 <- 50
  covariate <- c(rep(0, 10), seq_len(n - 10) / (n-10))
  true_variances0 <- 20 * covariate + 0.3
  # Based on the first equation of section 3 in the 2004 limma paper by Smyth
  true_variance <- 1/(rchisq(n = n, df = true_df0) / (true_df0 * true_variances0))
  observations <- t(vapply(seq_len(n), function(idx){
    rnorm(n = 5, mean = 0, sd = sqrt(true_variance[idx]))
  }, FUN.VALUE = rep(0.0, 5)))

  obs_var <- DelayedMatrixStats::rowVars(observations)
  df <- ncol(observations) - 1
  expect_silent(
    res_gp <- variance_prior(obs_var, df, covariate = covariate, abundance_trend = TRUE)
  )

  obs_var[33:40] <- 0
  expect_error(
    res_gp2 <- variance_prior(obs_var, df, covariate = covariate, abundance_trend = TRUE)
  )
})







