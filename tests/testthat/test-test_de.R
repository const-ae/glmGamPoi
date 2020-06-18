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


test_that("test_de works with contrast vector input", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol  =10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE),
                      cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- model.matrix(~ group + cont1 + cont2, data = annot)
  fit <- glm_gp(Y, design = design)
  res1 <- test_de(fit, contrast = cont2 - groupB)
  res2 <- test_de(fit, contrast = c(0, -1, 0, 1))
  expect_equal(res1, res2)

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
                  pseudobulk_by = sample,,
                  n_max = 4)

  res4 <- test_de(fit, reduced_design = ~ cont1 + 1,
                  full_design = ~ cont1 + cont2,
                  subset_to = cell_type == "Bcell",
                  pseudobulk_by = sample,
                  n_max = 4)
  res4
})


test_that("calling subset_to works", {

})

test_that("calling pseudo_bulk works", {

})


