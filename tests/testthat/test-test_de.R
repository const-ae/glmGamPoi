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

  test_de(fit, contrast = contr_mat)
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
    matrixStats::colMeans2(fit$model_matrix, rows = idx)
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


