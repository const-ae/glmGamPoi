test_that("test_de works", {
  set.seed(1)
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  reduced_des <- model.matrix(~ group + cont1, data = annot)
  fit <- glm_gp(Y, design = design)
  res <- test_de(fit, reduced_design = reduced_des)
  res2 <- test_de(fit, contrast = cont2)
  # Should be equal except for log-fold change column
  expect_equal(res[,-7], res2[,-7], tolerance = 1e-6)
  expect_equal(res$lfc, rep(NA, 30), tolerance = 1e-6)
  expect_equal(res2$lfc, fit$Beta[,"cont2"] / log(2), tolerance = 1e-6)

  design_wo_intercept <- model.matrix(~ group + cont1 + cont2 - 1, data = annot)
  red_design <- model.matrix( ~ cont1 + cont2 + 1, data = annot)
  fit_wo_intercept <- glm_gp(Y, design = design_wo_intercept)
  res3 <- test_de(fit_wo_intercept, reduced_design = red_design)
  res4 <- test_de(fit_wo_intercept, contrast = groupA - groupB)
  expect_equal(res3[,-7], res4[,-7], tolerance = 1e-6)

  res5 <- test_de(fit, contrast = "-groupB")
  expect_equal(res4, res5, tolerance = 0.1)
})

test_that("test_de works with 'fact' specifications", {
  set.seed(1)
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  reduced_des <- model.matrix(~ group + cont1, data = annot)
  fit <- glm_gp(Y, design = design)
  # Doesn't work because design is a matrix
  expect_error(test_de(fit, cond(group = "A") - cond(group = "B")))

  fit <- glm_gp(Y, design = ~ group + cont1 + cont2, col_data = annot)
  res1 <- test_de(fit, cond(group = "B") - cond(group = "A"))
  res2 <- test_de(fit, groupB)
  expect_equal(res1, res2)
})


test_that("test_de works with contrast vector input", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  fit <- glm_gp(Y, design = design)
  res1 <- test_de(fit, contrast = cont2 - groupB)
  res2 <- test_de(fit, contrast = c(0, -1, 0, 1))
  expect_equal(res1, res2)

  # Error on malformed contrast vector
  expect_error(test_de(fit, contrast = c(0, -1, 0)))

})

test_that("test_de works with a matrix as contrast", {
  Y <- matrix(rnbinom(n = 30 * 20, mu = 4, size = 0.3), nrow = 30, ncol  =20)
  annot <- data.frame(group = sample(c("A", "B", "C"), size = 20, replace = TRUE),
                      cont1 = rnorm(20), cont2 = rnorm(20, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2 -1, data = annot)
  fit <- glm_gp(Y, design = design)
  # res1 <- test_de(fit, contrast = cont2 - groupB)
  # res2 <- test_de(fit, contrast = c(0, -1, 0, 1))
  # expect_equal(res1, res2)
  contr_mat <- limma::makeContrasts(contrasts = c("groupA-groupB", "groupA - groupC", "groupB - groupC"), levels = colnames(fit$Beta))

  expect_no_error(test_de(fit, contrast = contr_mat))
})



test_that("test_de handles on_disk correctly", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  se <- SummarizedExperiment::SummarizedExperiment(Y_hdf5, colData = annot)
  # The actual check is difficult, because internally the res2 should be calculated on_disk
  fit1 <- glm_gp(se, design = ~ group + cont1 + cont2, on_disk = TRUE)
  res1 <- test_de(fit1, reduced_design = ~ 1)
  fit2 <- glm_gp(se, design = ~ group + cont1 + cont2, on_disk = FALSE)
  res2 <- test_de(fit2, reduced_design = ~ 1)
  expect_equal(res1, res2)

})



test_that("NSE works", {
  Y <- matrix(rnbinom(n = 30 * 100, mu = 4, size = 0.3), nrow = 30, ncol = 100)
  annot <- data.frame(sample = sample(c("A", "B", "C", "D"), size = 100, replace = TRUE),
                      cont1 = rnorm(100), cont2 = rnorm(100, mean = 30),
                      cell_type = sample(c("Tcell", "Bcell", "Makrophages"), size = 100, replace = TRUE))
  annot$condition <- ifelse(annot$sample %in% c("A", "B"), "ctrl", "treated")
  head(annot)
  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = annot)
  fit <- glm_gp(se, design = ~ condition + cont1 + cont2 + cell_type - 1, overdispersion = FALSE, size_factors = FALSE)

  expect_no_error({
    res <- test_de(fit, conditionctrl - conditiontreated)
    res2 <- test_de(fit, conditionctrl - conditiontreated,
                    full_design =  ~ condition + cont1 + cont2 - 1,
                    subset_to = cell_type == "Tcell",
                    n_max = 4)
    res3 <- test_de(fit, conditionctrl - conditiontreated,
                    full_design =  ~ condition + cont1 - 1,
                    subset_to = cell_type == "Tcell",
                    pseudobulk_by = sample,
                    n_max = 4)
    res4 <- test_de(fit, reduced_design = ~ cont1 + 1,
                    full_design = ~ cont1 + cont2,
                    subset_to = cell_type == "Bcell",
                    pseudobulk_by = sample,
                    n_max = 4)
  })
})


test_that("Pseudo bulk produces same results as making it manually", {
  Y <- matrix(rnbinom(n = 30 * 100, mu = 4, size = 0.3), nrow = 30, ncol = 100)
  annot <- data.frame(sample = sample(LETTERS[1:6], size = 100, replace = TRUE),
                      cont1 = rnorm(100), cont2 = rnorm(100, mean = 30),
                      cell_type = sample(c("Tcell", "Bcell", "Makrophages"), size = 100, replace = TRUE))
  annot$condition <- ifelse(annot$sample %in% c("A", "B", "C"), "ctrl", "treated")
  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = annot)

  fit <- glm_gp(se, design = ~ condition + cont1 + cont2 - 1, ridge_penalty = 5)
  res <- test_de(fit, reduced_design = ~ cont1 + cont2 + 1,
                 pseudobulk_by = sample)

  splitter <- split(seq_len(ncol(se)), SummarizedExperiment::colData(se)$sample)
  pseudobulk_mat <- do.call(cbind, lapply(splitter, function(idx){
    matrixStats::rowSums2(assay(se), cols = idx)
  }))
  pseudo_design_mat <- do.call(rbind, lapply(splitter, function(idx){
    matrixStats::colMeans2(fit$model_matrix, rows = idx, useNames = FALSE)
  }))


  fit2 <- glm_gp(pseudobulk_mat, design = pseudo_design_mat, ridge_penalty = 5)
  res2 <- test_de(fit2, contrast = Coef_1 - Coef_2)

  # Equal except for lfc column because of contrast vs reduced_design stuff:
  expect_equal(res[,-7], res2[,-7])

  res3 <- test_pseudobulk(se, design = ~ condition + cont1 + cont2 - 1,
                          aggregate_cells_by = sample,
                          contrast = conditionctrl - conditiontreated,
                          ridge_penalty = 5)
  expect_equal(res2, res3)

  res4 <- test_pseudobulk(se, design = ~ condition + cont1 + cont2 - 1,
                          aggregate_cells_by = sample,
                          contrast = conditionctrl - conditiontreated)

  res5 <- test_pseudobulk(se, design = ~ condition + cont1 + cont2,
                          reference_level = "treated",
                          aggregate_cells_by = sample,
                          contrast = conditionctrl)
  # They are as good as identical
  lapply(c("pval", "adj_pval", "f_statistic", "lfc"), function(col){
    expect_gt(cor(res4[[col]], res5[[col]]), 0.99)
  })

})


test_that("pseudobulk works with a matrix as contrast", {
  Y <- matrix(rnbinom(n = 30 * 100, mu = 4, size = 0.3), nrow = 30, ncol = 100)
  annot <- data.frame(sample = sample(LETTERS[1:6], size = 100, replace = TRUE),
                      cont1 = rnorm(100), cont2 = rnorm(100, mean = 30),
                      cell_type = sample(c("Tcell", "Bcell", "Makrophages"), size = 100, replace = TRUE))
  annot$condition <- ifelse(annot$sample %in% c("A", "B", "C"), "ctrl", "treated")
  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = annot)
  fit <- glm_gp(se, design = ~ condition + cont1 + cont2 - 1, ridge_penalty = 5)
  contr_mat <- limma::makeContrasts(contrasts = c("conditiontreated - conditionctrl", "cont1 - conditionctrl", "cont1 - cont2"), levels = colnames(fit$Beta))

  res <- test_de(fit, contrast = contr_mat,
                 pseudobulk_by = sample)

  splitter <- split(seq_len(ncol(se)), SummarizedExperiment::colData(se)$sample)
  pseudobulk_mat <- do.call(cbind, lapply(splitter, function(idx){
    matrixStats::rowSums2(assay(se), cols = idx)
  }))
  pseudo_design_mat <- do.call(rbind, lapply(splitter, function(idx){
    matrixStats::colMeans2(fit$model_matrix, rows = idx)
  }))


  fit2 <- glm_gp(pseudobulk_mat, design = pseudo_design_mat, ridge_penalty = 5)
  res2 <- test_de(fit2, contrast = contr_mat)

  expect_equal(res, res2)
})



test_that("offset is correctly propagated in test_de()", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  reduced_des <- model.matrix(~ group + cont1, data = annot)

  fit1 <- glm_gp(Y, design = design)
  offset <- matrix(log(fit1$size_factors), byrow = TRUE, nrow = nrow(Y), ncol = ncol(Y))
  fit2 <- glm_gp(Y, design = design, offset = offset, size_factors = FALSE)

  expect_equal(fit1[-6], fit2[-6])

  res1 <- test_de(fit1, reduced_design = reduced_des)
  res2 <- test_de(fit2, reduced_design = reduced_des)

  expect_equal(res1, res2)
})


test_that("likelihood ratio test works with ridge_penalty", {

  set.seed(1)
  size <- 50
  group <- sample(LETTERS[1:4], size = size, replace = TRUE)
  mu <- c(A = 2, B = 4, C = 10, D = 2.5)
  y <- rpois(n = size, lambda = mu[group])
  fit1 <- glm_gp(y, group, ridge_penalty = 0, overdispersion = 0)
  res1 <- test_de(fit1, B - D)

  fit2 <- glm_gp(y, group, ridge_penalty = 1e-6, overdispersion = 0)
  res2 <- test_de(fit2, B - D)

  expect_equal(fit1[- which(names(fit1) == "ridge_penalty")],
               fit2[- which(names(fit2) == "ridge_penalty")])
  expect_equal(res1, res2)

  fit3 <- glm_gp(y, group, ridge_penalty = 5, overdispersion = 0)
  res3 <- test_de(fit3, B - D)

  expect_true(all(fit3$Beta < fit1$Beta))
  expect_gt(res3$f_statistic, 0)
  expect_gt(res3$pval, res2$pval)


  # Should be similar but not identical
  fit3 <- glm_gp(y, group, ridge_penalty = 10)
  fit4 <- glm_gp(y, group, ridge_penalty = 10, reference_level = "A")
  fit3$Beta
  fit4$Beta
  test_de(fit3, B - D)
  test_de(fit4, B_vs_A - D_vs_A)


})


test_that("ridge_penalty is adapted if reduced_design is specified", {

  size <- 100
  y <- rnbinom(n = size, mu = 5, size = 1 / 3.1)
  df <- data.frame(group = sample(LETTERS[1:2], size, replace = TRUE),
                   cont = rnorm(size))

  fit <- glm_gp(y, design = ~ group + cont + group:cont, col_data = df)
  res1 <- test_de(fit, contrast = `groupB:cont`)
  res2 <- test_de(fit, reduced_design = ~ group + cont)
  expect_equal(res1[,-which(names(res1) == "lfc")],
               res2[,-which(names(res2) == "lfc")],
               tolerance = 1e-7)

  fit_ridge <- glm_gp(y, design = ~ group + cont + group:cont, col_data = df,
                      ridge_penalty = 1.2, reference_level = "B")
  res1 <- test_de(fit_ridge, contrast = `groupA:cont`)
  res2 <- test_de(fit_ridge, reduced_design = ~ group + cont)
  expect_equal(res1[,-which(names(res1) == "lfc")],
               res2[,-which(names(res2) == "lfc")],
               tolerance = 1e-6)


})


test_that("ignore_degeneracy is propagated", {
  y <- rnbinom(n = 30, mu = 7, size = 1/2.2)
  X <- cbind(matrix(rnorm(n = 30 * 3), nrow = 30, ncol = 3), 0)
  Xmod <- X
  attr(Xmod, "ignore_degeneracy") <- TRUE

  expect_error(glm_gp(y, design = X, ridge_penalty = 0.1))

  # These succeed
  fit <- glm_gp(y, design = Xmod, ridge_penalty = 0.1)
  test_de(fit, contrast = c(0, 1, 0, 0))
})


test_that("test_de works without shrinkage", {
  set.seed(1)
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  fit <- glm_gp(Y, design = ~ group + cont1, col_data = annot)
  res <- test_de(fit, reduced_design = ~ group)

  fit2 <- glm_gp(Y, design = ~ group + cont1, col_data = annot, overdispersion_shrinkage = FALSE)
  res2 <- test_de(fit2, reduced_design = ~ group)

  expect_equal(colnames(res), colnames(res2))
  expect_equal(res$name, res2$name)
  expect_equal(res$df1, res2$df1)
  expect_equal(res2$df2, rep(NA_real_, 30))
  expect_gt(cor(res$pval, res2$pval), 0.9)
})


test_that("test_de works with compute_lfc_se", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 3000, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  cntrst <- c(0, -1, 0, 0)
  fit <- glm_gp(Y, design = design, overdispersion_shrinkage = TRUE, size_factors = "ratio")
  res <- test_de(fit, contrast = cntrst,
                 compute_lfc_se=TRUE)
  pred <- predict(fit, se.fit=TRUE, newdata=matrix(cntrst, nrow=1))

  expect_equal(res$lfc_se, pred$se.fit[,1] / log(2))

  res2 <- test_de(fit, reduced_design = model.matrix(~ cont1 + cont2, data = annot),
                 compute_lfc_se=TRUE)

  expect_equal(res$pval, res2$pval)
  expect_true(all(is.na(res2$lfc)))
  expect_true(all(is.na(res2$lfc_se)))

  res3 <- test_de(fit, contrast = cntrst)
  res4 <- test_de(fit, reduced_design = model.matrix(~ cont1 + cont2, data = annot))
  expect_equal(res$pval, res3$pval)
  expect_false("lfc_se" %in% colnames(res3))
  expect_equal(res$lfc, res3$lfc)

  expect_equal(res$pval, res4$pval)
  expect_true(all(is.na(res4$lfc)))
  expect_false("lfc_se" %in% colnames(res3))
})

test_that("test_de works with degenerate designs", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 3000, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group, data = annot)
  design <- cbind(design, design[,2])
  attr(design, "ignore_degeneracy") <- TRUE
  cntrst <- matrix(c(0, 1, 0), nrow = 1)
  fit <- glm_gp(Y, design = design, overdispersion_shrinkage = TRUE, size_factors = "ratio",
                ridge_penalty = 1e-3)
  skip("test_de with compute_lfc_se cannot handle degenerate designs")
  res <- test_de(fit, contrast = c(cntrst), compute_lfc_se=TRUE)
  expect_true(all(res$lfc_se > 10))

  # res$pval_alt <- with(res, 2 * pnorm(abs(lfc / lfc_se), lower.tail = FALSE))
  # head(res)
  # expect_equal(0.3 * fit$Beta[,2], fit$Beta[,3], tolerance = 1e-4)
  # pred <- predict(fit, se.fit=TRUE, newdata=matrix(cntrst, nrow=1))
  #
  # # Wald Test P-value
  # plot(res$pval, res$pval_alt); abline(0,1)
  #
  # gene_idx <- 1
  # disp <- fit$overdispersions[gene_idx]
  # mu <- fit$Mu[gene_idx,]
  # weights <- mu / (1 + mu * disp)
  # weighted_Design <-  fit$model_matrix * sqrt(weights)
  # ridge_penalty <- diag(1e-3, nrow = ncol(fit$model_matrix))
  #
  # # Below is an optimized formula to calculate the variance covariance matrix (X'X)^-1
  # # See https://stats.stackexchange.com/a/407734/130486
  # Rinv <- qr.solve(qr.R(qr(weighted_Design)))
  # lhs <- cntrst %*% Rinv
  # sqrt(matrixStats::rowSums2(lhs^2))
  #
  # # This is the explicit formula to calculate the variance covariance matrix with ridge.
  # # Formula is from "Wald test" section page 18 of  DESeq2 paper (Love et al., 2014)
  # XtwX <- t(weighted_Design) %*% weighted_Design
  # XtwX_ridge_inv <- solve(XtwX + nrow(weighted_Design) * t(ridge_penalty) %*% ridge_penalty)
  # cov_mat <- XtwX_ridge_inv %*% XtwX %*% XtwX_ridge_inv
  # sqrt(diag(cntrst %*% cov_mat %*% t(cntrst)))
  #
  # # This is an optimized implementation that is numerically more robust
  # Xwave <- rbind(weighted_Design, sqrt(nrow(weighted_Design)) * ridge_penalty)
  # Rinv <- qr.solve(qr.R(qr(Xwave)))
  # lhs <- cntrst %*% (Rinv %*% t(Rinv)) %*% t(weighted_Design)
  # sqrt(matrixStats::rowSums2(lhs^2))
  #
  # # Fourth option
  # Qr <- qr(Xwave)
  # p1 <- seq_len(Qr$rank)
  # cov_mat2 <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  # sqrt(diag(cntrst %*% cov_mat2 %*% t(cntrst)))
  #
  # lm_fit <- treelabel:::modern_glm(1:10, design = design, family = gaussian())
  # lm_fit$coefficients
  # treelabel:::test_diff(lm_fit, cntrst)
  #
  # ################# Linear model
  # design <- model.matrix(~ group, data = annot)
  # design <- cbind(design, design[,2])
  # cntrst <- matrix(c(0, 1, 0), nrow = 1)
  # ridge_penalty <- diag(1e-3, nrow = ncol(design))
  # y <- rnorm(10)
  # lm_fit <- lm(y ~ design - 1)
  # coef(lm_fit)
  # c(solve(t(design) %*% design + ridge_penalty) %*% t(design) %*% y)
  # lm_summary <- summary(lm(y ~ design - 1))
  # lm_summary$sigma^2 * lm_summary$cov.unscaled
  # lm_summary$sigma^2 * solve(t(design) %*% design + ridge_penalty)
  # cov_mat <- solve(t(design) %*% design + ridge_penalty) %*% t(design) %*% diag(residuals(lm_fit)^2) %*% design %*% solve(t(design) %*% design + ridge_penalty)
  # cntrst %*% cov_mat %*% t(cntrst)
  # sandwich::vcovHC(lm_fit, type = "HC0")
  #
  # V <- svd(design)$v
  # mod_design <- design %*% V
  # cov_mat2 <- solve(t(mod_design) %*% mod_design + ridge_penalty) %*% t(mod_design) %*% diag(residuals(lm_fit)^2) %*% mod_design %*% solve(t(mod_design) %*% mod_design + ridge_penalty)
  # V %*% cov_mat2 %*% t(V)
})


# test_that("glmGamPoi results are similar to DESeq2",{
#   Y <- matrix(rnbinom(n = 30 * 10, mu = 3000, size = 0.3), nrow = 30, ncol  =10)
#   annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
#                       cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
#   design <- model.matrix(~ group, data = annot)
#   cntrst <- c(0, -1)
#   fit <- glm_gp(Y, design = design, overdispersion_shrinkage = TRUE, size_factors = "ratio")
#   res <- test_de(fit, contrast = cntrst,
#                  compute_lfc_se=TRUE)
#   pred <- predict(fit, se.fit=TRUE, newdata=matrix(cntrst, nrow=1))
#
#   pval <- pnorm(res$lfc / res$lfc_se)
#   pval <- 2 * pmin(pval, 1-pval)
#   cor(res$pval, pval)
#   plot(res$pval, pval, log = "xy"); abline(0,1)
#
#   # # Compare against DESeq2
#   dds <- DESeq2::DESeqDataSetFromMatrix(Y, annot, design = ~ group)
#   # dds <- DESeq2::DESeq(dds)
#   dds <- DESeq2::estimateSizeFactors(dds)
#   expect_equal(cor(fit$size_factors, dds$sizeFactor), 1)
#
#   # dds <- DESeq2:::estimateDispersions.DESeqDataSet(dds, fitType = "glmGamPoi")
#   dds <- DESeq2::estimateDispersionsGeneEst(dds, type = "glmGamPoi", niter = 2)
#   expect_equal(fit$overdispersions, rowData(dds)$dispGeneEst, tolerance = 1e-4)
#   dds <- DESeq2::estimateDispersionsFit(dds, fitType = "glmGamPoi")
#   rowData(dds)$dispFit <- fit$overdispersion_shrinkage_list$dispersion_trend
#   dds <- DESeq2::estimateDispersionsMAP(dds, type = "glmGamPoi", dispPriorVar = 1e6)
#   dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
#   res_deseq <- DESeq2::results(dds, contrast = cntrst, minmu = 1e-6)
#   plot(res$lfc_se, res_deseq$lfcSE); abline(0,1)
#
#
#   pval <- 2 * pnorm(abs(res_deseq$log2FoldChange/ res_deseq$lfcSE), lower.tail = FALSE)
#   plot(pval, res_deseq$pvalue); abline(0,1)
#   pval2 <- 2 * pnorm(abs(res$lfc/ res$lfc_se), lower.tail = FALSE)
#   plot(pval2, res_deseq$pvalue, col = as.factor(abs(res_deseq$lfcSE - res$lfc_se) < 1e-4)); abline(0,1)
#
#
#   # Test with custom ridge penalty. DESeq doesn't support arbitrary lambda values, so
#   # patch by manually changing the values through the debugger.
#   lambda_global <- c(0.2, 3)
#   lambda_deseq <- lambda_global^2 * ncol(Y)
#   debugonce(DESeq2:::getContrast)
#   res_deseq2 <- DESeq2::results(dds, contrast = cntrst, minmu = 1e-6)
#   fit2 <- glm_gp(Y, design = design, overdispersion_shrinkage = TRUE, size_factors = "ratio",
#                  ridge_penalty = lambda_global)
#   res2 <- test_de(fit2, contrast = cntrst, compute_lfc_se=TRUE)
#   # Results agree well enough
#   plot(res2$lfc_se, res_deseq2$lfcSE); abline(0,1)
#
#
#   # # Compare against edgeR
#   edgeR_data <- edgeR::DGEList(Y)
#   edgeR_data <- edgeR::calcNormFactors(edgeR_data, method =  "RLE")
#   edgeR_data <- edgeR::estimateDisp(edgeR_data, design)
#   edgeR_fit <- edgeR::glmFit(edgeR_data, design = design)
#   edgeR_fit <- edgeR::glmLRT(edgeR_fit, contrast = cntrst)
#   res_edgeR <- edgeR::topTags(edgeR_fit, sort.by = "none", n = nrow(Y))
#   plot(res$lfc, res_edgeR$table$logFC); abline(0,1)
#   plot(res$pval, res_edgeR$table$PValue); abline(0,1)
# })
