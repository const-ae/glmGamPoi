

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


test_that("variance prior estimation works", {
  n <- 1000
  true_df0 <- 5
  true_variance0 <- 1.3
  # Based on the first equation of section 3 in the 2004 limma paper by Smyth
  true_variance <- 1/(rchisq(n = n, df = true_df0) / (true_df0 * true_variance0))
  observations <- t(vapply(seq_len(n), function(idx){
    rnorm(n = 5, mean = 0, sd = sqrt(true_variance[idx]))
  }, FUN.VALUE = rep(0.0, 5)))

  obs_var <- matrixStats::rowVars(observations)
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

  obs_var <- matrixStats::rowVars(observations)
  df <- ncol(observations) - 1
  res_limma <- limma:::squeezeVar(obs_var, df = df, covariate = covariate)
  res_gp <- variance_prior(obs_var, df, covariate = covariate)
  res_gp2 <- variance_prior(obs_var, df, covariate = NULL)


  loglikelihood_gp <- sum(df(obs_var/res_gp$variance0, df1=df, df2=res_gp$df0, log=TRUE) - log(res_gp$variance0))
  loglikelihood_limma <- sum(df(obs_var/res_limma$var.prior, df1=df, df2=res_limma$df.prior, log=TRUE) - log(res_limma$var.prior))
  expect_gt(loglikelihood_gp, loglikelihood_limma)

})



test_that("shrinkage makes sense", {

  se <- HighlyReplicatedRNASeq::Schurch16()
  colData(se)$cont <- rnorm(n = ncol(se))
  fit <- glm_gp(se, design = ~ condition + cont, verbose = TRUE)

  res <- gampoi_test_qlr(se, fit, reduced_design = ~ cont, verbose = TRUE)

  hist(res$pval)

  gene_means <- rowMeans(fit$Mu)

  plot(gene_means, fit$overdispersions, log = "xy", cex = 0.7)
  # points(gene_means, fit$overdispersion_shrunken, col = "blue", pch = 16, cex = 0.4)
  # lines(sort(gene_means), fit$overdispersion_trend[order(gene_means)], col = "red")

  # s2 <- rowSums(residuals.glmGamPoi(fit, assay(se))^2) / (ncol(se) - ncol(fit$model_matrix))

  plot(gene_means, fit$ql_disp_estimate^(1/4), log = "x")
  lines(sort(gene_means), fit$ql_disp_trend[order(gene_means)]^(1/4), col ="red")
  points(gene_means, fit$ql_disp_shrunk^(1/4), col ="blue", pch  =16 , cex = 0.3)

})













