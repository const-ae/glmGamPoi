test_that("glm_gp works for simple cases", {
  y <- rpois(n = 10, lambda = 3)
  expect_silent({
    tmp <- glm_gp(y, design = ~ 1)
  })
  tmp <- glm_gp(y, design = ~ 1, size_factors = FALSE)
  expect_equal(tmp$size_factors, rep(1, times = 10))
  expect_equal(c(tmp$model_matrix), rep(1, times = 10))
  expect_equal(c(tmp$Mu), rep(mean(y), times = 10))
  expect_equal(c(tmp$Beta), log(mean(y)))
})

test_that("glm_gp() can be called like glm(): formula first", {

  y <- c(1, 3, 6)
  expect_equal(glm_gp(y ~ 1), glm_gp(y))

  Y <- matrix(rpois(n = 4 * 9, lambda = 2), nrow = 4, ncol = 9)
  expect_equal(glm_gp(Y ~ 1), glm_gp(Y))

  a <- rnorm(9)
  df <- data.frame(b = runif(9),
                   y = rpois(9, lambda = 3))
  y <- df$y + 10
  expect_equal(glm_gp(Y ~ a + b, col_data = df), glm_gp(Y, ~ a + b, col_data = df))
  expect_error(expect_equal(glm_gp(y ~ a + b, col_data = df), glm_gp(y, ~ a + b, col_data = df)))
  expect_equal(glm_gp(y ~ a + b, col_data = df), glm_gp(df$y, ~ a + b, col_data = df))

  expect_error(glm_gp(~ a + b))
  expect_error(glm_gp(y ~ a,  ~ a))
  glm_gp(y ~ ., col_data = df)$design_formula

})


test_that("glm_gp can handle complicated formulas", {

  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(100))
  colnames(Y) <- paste0("person_", seq_len(50))

  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
  expect_equal(colnames(fit$Beta), colnames(fit$model_matrix))
  expect_equal(rownames(fit$Beta), rownames(Y))
  expect_equal(colnames(fit$Mu), colnames(Y))
  expect_equal(rownames(fit$Mu), rownames(Y))
  expect_equal(rownames(fit$model_matrix), colnames(Y))
  expect_equal(names(fit$overdispersions), rownames(Y))
  expect_equal(names(fit$size_factors), colnames(Y))
})

test_that("glm_gp can handle overdispersion = 'global'", {
  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(100))
  colnames(Y) <- paste0("person_", seq_len(50))

  fit <- glm_gp(Y, overdispersion = "global")
  res <- overdispersion_mle(Y, fit$Mu, global_estimate = TRUE)
  expect_equal(unname(fit$overdispersions[1]), res$estimate)
  expect_equal(unname(fit$overdispersion_shrinkage_list$dispersion_trend), unname(fit$overdispersions))
  expect_equal(unname(fit$overdispersion_shrinkage_list$ql_disp_estimate[1]), 1)
})


test_that("handle_design_parameter throws error for degenerate design matrix", {

  df <- data.frame(group = sample(LETTERS[1:3], size = 20, replace = TRUE),
             cont = rnorm(20))
  df$group_copy <- df$group
  data <- matrix(0, ncol = 20, nrow = 1)
  expect_error(handle_design_parameter(~ group + group_copy + cont, data, col_data = df, reference_level = NULL))

  model_matrix <- model.matrix(~ group + cont, data = df)
  model_matrix <- cbind(model_matrix, 0)
  expect_error(handle_design_parameter(model_matrix, data, col_data = df, reference_level = NULL))

})

test_that("handle_design_parameter throws error for design matrix with NA's", {
  data <- matrix(0, ncol = 20, nrow = 1)
  df <- data.frame(var=c(1:19, NA))
  expect_error(handle_design_parameter(~ var, data, col_data = df, reference_level = NULL))
})

test_that("glm_gp can handle no-row input", {

  Y <- matrix(numeric(0), nrow = 0, ncol = 10)
  tmp <- glm_gp(Y, size_factors = FALSE)
  expect_equal(dim(tmp$Beta), c(0, 1))
  expect_equal(dim(tmp$Mu), c(0, 10))
  expect_equal(tmp$size_factors, rep(1, times = 10))
  expect_equal(c(tmp$model_matrix), rep(1, times = 10))


})

test_that("glm_gp throws error on no-col input", {

  Y <- matrix(numeric(0), nrow = 3, ncol = 0)
  expect_error(glm_gp(Y, size_factors = FALSE))


})


test_that("glm_gp produces appropriate error message for NA input", {

  y <- c(NA, rpois(n = 9, lambda = 3), NA, NA, NA, NA, NA, NA)
  expect_error(
    tmp <- glm_gp(y, size_factors = FALSE)
  )
})

test_that("glm_gp produces appropriate error message for Infinite input", {

  y <- c(-Inf, rpois(n = 9, lambda = 3))
  expect_error(
    tmp <- glm_gp(y, size_factors = FALSE)
  )

})

test_that("glm_gp produces appropriate error message for negative input", {

  y <- c(-3, rpois(n = 9, lambda = 3))
  expect_error(
    tmp <- glm_gp(y, size_factors = FALSE)
  )

})

test_that("glm_gp works for subsample = Inf", {
  y <- rpois(n = 10, lambda = 3)
  tmp <- glm_gp(y, design = ~ 1, size_factors = FALSE)
  tmpInf <- glm_gp(y, design = ~ 1, size_factors = FALSE, subsample = Inf)
  expect_equal(tmp, tmpInf)
})


test_that("glm_gp can handle on_disk parameter", {

  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))
  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)
  dimnames(Y_hdf5) <- dimnames(Y)

  fit_in_memory <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data, on_disk = FALSE)
  fit_on_disk <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data, on_disk = TRUE)

  expect_equal(fit_in_memory[!names(fit_in_memory) %in% c("Mu", "data", "Offset")], fit_on_disk[!names(fit_on_disk) %in% c("Mu", "data", "Offset")])
  expect_equal(c(fit_in_memory$Mu), c(fit_on_disk$Mu))
  expect_s4_class(fit_on_disk$Mu, "DelayedArray")
  expect_equal(c(assay(fit_in_memory$data)), c(assay(fit_on_disk$data)))
  expect_equal(class(fit_in_memory$data), class(fit_on_disk$data))
  expect_s4_class(assay(fit_on_disk$data), "DelayedArray")
  expect_equal(c(fit_in_memory$Offset), c(fit_on_disk$Offset))
  expect_s4_class(fit_on_disk$Offset, "DelayedArray")


  fit_on_disk2 <- glm_gp(Y_hdf5, design = ~ fav_food + city + age, col_data = data)
  expect_equal(fit_on_disk[!names(fit_in_memory) %in% c("Mu", "data", "Offset")], fit_on_disk2[!names(fit_on_disk2) %in% c("Mu", "data", "Offset")])
  expect_s4_class(fit_on_disk2$Mu, "DelayedArray")
  expect_s4_class(assay(fit_on_disk2$data), "DelayedArray")
  expect_s4_class(fit_on_disk2$Offset, "DelayedArray")

  fit_in_memory2 <- glm_gp(Y_hdf5, design = ~ fav_food + city + age, col_data = data, on_disk = FALSE)
  expect_equal(fit_in_memory, fit_in_memory2)

})




test_that("glm_gp can handle SummarizedExperiment correctly", {


  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))

  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = data)

  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
  fit_se <- glm_gp(se, design = ~ fav_food + city + age)

  expect_equal(fit, fit_se)
})


test_that("glm_gp can handle sparse input correctly", {


  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))
  Y_sp <- as(Y, "dgCMatrix")

  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
  expect_error({
    glm_gp(Y_sp, design = ~ fav_food + city + age, col_data = data)
  })
  fit_sp1 <- glm_gp(Y_sp, design = ~ fav_food + city + age, col_data = data, on_disk = FALSE)
  fit_sp2 <- glm_gp(Y_sp, design = ~ fav_food + city + age, col_data = data, on_disk = TRUE)

  expect_equal(fit, fit_sp1)

  expect_equal(fit[!names(fit) %in% c("Mu", "data", "Offset")], fit_sp2[!names(fit_sp2) %in% c("Mu", "data", "Offset")])
  expect_equal(c(fit$Mu), c(fit_sp2$Mu))
  expect_s4_class(fit_sp2$Mu, "DelayedArray")
  expect_equal(c(assay(fit$data)), c(assay(fit_sp2$data)))
  expect_equal(class(fit$data), class(fit_sp2$data))
  expect_s4_class(assay(fit_sp2$data), "DelayedArray")
  expect_equal(c(fit$Offset), c(fit_sp2$Offset))
  expect_s4_class(fit_sp2$Offset, "DelayedArray")
})

test_that("glm_gp can handle intercept model", {


  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 150, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 150, replace = TRUE),
                     age = rnorm(n = 150, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 150, mu = 3, size = 1/3.1), nrow = 10, ncol = 150)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(150))

  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = data)

  expect_silent({
    fit1 <- glm_gp(se, design = ~ fav_food + city + age)
    fit2 <- glm_gp(Y, design = ~ 1)
    fit3 <- glm_gp(se, design = ~ fav_food + city + age, overdispersion = 0)
  })


})


test_that("glm_gp can handle design parameter of type vector", {

  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))

  fit_intercept <- glm_gp(Y, design = ~ 1)
  expect_error(glm_gp(Y, design = 1))

  sample_assignment <- rep(1, times = 50)
  fit_vec <- glm_gp(Y, design = sample_assignment)

  expect_equal(c(fit_intercept$Beta), c(fit_vec$Beta))
  expect_equal(colnames(fit_vec$Beta), "1")

  sample_assignment <- sample(c("a", "b", "c"), size = 50, replace = TRUE)
  design_mat <- model.matrix(~ x_ - 1, data = data.frame(x_ = sample_assignment))
  colnames(design_mat) <- c("a", "b", "c")
  fit_vec <- glm_gp(Y, design = sample_assignment)
  fit_mat <- glm_gp(Y, design = design_mat)
  # Ignore design_formula
  expect_equal(fit_vec[-c(10, 11)], fit_mat[-10])


  sample_assignment <- sample(c("a", "b", "c"), size = 50, replace = TRUE)
  design_mat <- model.matrix(~ x_, data = data.frame(x_ = relevel(as.factor(sample_assignment), ref = "b")))
  fit_vec <- glm_gp(Y, design = sample_assignment, reference_level = "b")
  expect_error(glm_gp(Y, design = design_mat, reference_level = "b"))
  fit_mat <- glm_gp(Y, design = design_mat)
  fit_formula <- glm_gp(Y, design = ~ x_, col_data = data.frame(x_ = sample_assignment),
                        reference_level = "b")
  expect_equal(unname(fit_vec$Beta), unname(fit_mat$Beta))
  expect_equal(colnames(fit_vec$Beta), c("Intercept", "a_vs_b", "c_vs_b"))
  expect_equal(unname(fit_vec$Beta), unname(fit_formula$Beta))
  expect_equal(colnames(fit_formula$Beta), c("Intercept", "x_a", "x_c"))
})



test_that("glm_gp can handle design formula correctly", {

  coldata <- data.frame(condition = c(rep("A", 4), rep("B", 3)), stringsAsFactors = FALSE)
  Y <- matrix(numeric(0), ncol = 7, nrow = 0)
  fit <- glm_gp(Y, design = ~ condition, col_data = coldata, size_factors = FALSE)
  expect_equal(c(fit$model_matrix), c(model.matrix(~ condition, coldata)))

})


test_that("glm_gp gives error for too large col_data ", {

  coldata <- data.frame(condition = c(rep("A", 4), rep("B", 4)), stringsAsFactors = FALSE)
  rownames(coldata) <- paste0("sample_", sample(1:8))
  Y <- matrix(numeric(0), ncol = 7, nrow = 0)
  colnames(Y) <- paste0("sample_", 1:7)
  expect_error(glm_gp(Y, design = ~ condition, col_data = coldata, size_factors = FALSE))

})

test_that("glm_gp gives error for wrong input data ", {
  expect_error(glm_gp(data = c("hello", "world")))
  sp_mat <- as(matrix(1:10, nrow = 1), "dgCMatrix")
  expect_error(glm_gp(data = sp_mat))
})

test_that("glm_gp can handle mismatching assay vs model matrix names ", {

  coldata <- data.frame(condition = c(rep("A", 4), rep("B", 3)), stringsAsFactors = FALSE)
  rownames(coldata) <- paste0("sample_", sample(1:7))
  Y <- matrix(numeric(0), ncol = 7, nrow = 0)
  colnames(Y) <- paste0("sample_", 1:7)
  # Fails because SummarizedExperiment now checks that data and col_data match
  # expect_error(glm_gp(Y, design = ~ condition, col_data = coldata, size_factors = FALSE))

  # Works because this is only reordering the model matrix and doesn't have a model matrix
  mm <- model.matrix(~ condition, coldata)
  fit <- glm_gp(Y, design = mm)
  expect_equal(colnames(fit$data), paste0("sample_", 1:7))
  expect_equal(rownames(fit$model_matrix), paste0("sample_", 1:7))
})


test_that("glm_gp doesn't copy the data", {

  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))
  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)
  dimnames(Y_hdf5) <- dimnames(Y)

  fit_on_disk <- glm_gp(Y_hdf5, design = ~ fav_food + city + age, col_data = data)
  expect_s4_class(assay(fit_on_disk$data), "DelayedArray")
  expect_equal(assay(fit_on_disk$data)@seed@seed@filepath, Y_hdf5@seed@seed@filepath)


  # fit_in_memory <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data, on_disk = FALSE)
  # # If you ignore the wrapper for the second, the two objects are identical
  # .Internal(inspect(Y))
  # .Internal(inspect(assay(fit_in_memory$data)))
})



test_that("lte_n_equal_rows and get_row_groups works", {
  df <- data.frame(v1 = rep(LETTERS[1:3], each = 10),
                   v2 = rep(LETTERS[1:3], each = 10),
                   cont = rnorm(30), stringsAsFactors = TRUE)
  X1 <- model.matrix(~ v1 - 1, data = df)
  X2 <- model.matrix(~ v1 + 1, data = df)
  X3 <- model.matrix(~ v1 + v2 + 1, data = df)
  X4 <- model.matrix(~ v1 + v2 + cont, data = df)

  expect_true(lte_n_equal_rows(X1, n = ncol(X1)))
  expect_false(lte_n_equal_rows(X1, n = 2))
  expect_true(lte_n_equal_rows(X1, n = 4))
  expect_equal(get_row_groups(X1, n_groups = ncol(X1)),
               as.numeric(df$v1))
  expect_equal(get_groups_for_model_matrix(X1),
               as.numeric(df$v1))

  expect_true(lte_n_equal_rows(X2, n = 3))
  expect_false(lte_n_equal_rows(X2, n = 2))
  expect_equal(get_row_groups(X2, n_groups = ncol(X2)),
               as.numeric(df$v1))
  expect_equal(get_groups_for_model_matrix(X2),
               as.numeric(df$v1))


  expect_true(lte_n_equal_rows(X3, n = ncol(X3)))
  expect_true(lte_n_equal_rows(X3, n = 3))
  expect_false(lte_n_equal_rows(X3, n = 2))
  expect_equal(get_row_groups(X3, n_groups = ncol(X3)),
               as.numeric(df$v1))
  expect_equal(get_groups_for_model_matrix(X3),
               as.numeric(df$v1))


  expect_false(lte_n_equal_rows(X4, n = 20))
  expect_equal(get_row_groups(X4, n_groups = nrow(X4)), seq_len(nrow(df)))
  expect_null(get_groups_for_model_matrix(X4))

  # df2 <- data.frame(v1 = rep(LETTERS[1:3], each = 1000),
  #                  v2 = rep(LETTERS[1:3], each = 1000),
  #                  cont = rnorm(30))
  # X_large <- model.matrix(~ v1 + v2 + cont - 1, data = df2)
  #
  # bench::mark(
  #   lte_n_equal_rows(X_large, n = 10),
  #   DESeq2:::nOrMoreInCell(X_large, 10),
  #   check = FALSE
  # )

})



test_that("NA's produced by fitBeta_fisher_scoring don't cause problems", {

  skip_if_not(isNamespaceLoaded("devtools"))

  n_genes <- 100
  n_cells <- 500
  sf <- rchisq(n = n_cells, df = 5)
  sf <- sf / mean(sf)
  gene_means <- 10^runif(n = n_genes, min = -3, max = 3)

  Mu <- gene_means %*% t(sf)
  Y <- matrix(rnbinom(n = n_genes * n_cells, size = 1 / 0.1, mu = Mu), nrow = n_genes, ncol = n_cells)
  rownames(Y) <- paste0("Gene_", seq_len(n_genes))
  colnames(Y) <- paste("Cell_", seq_len(n_cells))
  cont <- rnorm(n_cells)
  cont2 <- rnorm(n_cells)

  fit_orig <- glm_gp(Y, design = ~ cont + cont2, size_factors = sf)


  # mockr::with_mock and mockery::stub don't properly work yet
  testthat::with_mock(
    # Simulated failure to converge
    fitBeta_fisher_scoring =  function(Y, model_matrix, exp_offset_matrix, thetas, beta_matSEXP, ridge_penalty_nl, tolerance, max_rel_mu_change, max_iter) {
      res <- .Call(`_glmGamPoi_fitBeta_fisher_scoring`, Y, model_matrix, exp_offset_matrix, thetas, beta_matSEXP, ridge_penalty_nl, tolerance, max_rel_mu_change, max_iter)
      if(! any(res$iter == 0)){
        res$beta_mat[c(5, 17), ] <- NA
        res$iter[c(5, 17)] <- 1000
        res$deviance[c(5, 17)] <- NaN
      }
      res
    },
    {
      # optim recovers problematic result
      # although they are not 100% identical
      expect_silent(fit <- glm_gp(Y, design = ~ cont + cont2, size_factors = sf))
      expect_equal(fit_orig$Beta[-c(5, 17),], fit$Beta[-c(5, 17),])
      expect_equal(fit_orig$Beta[c(5, 17),], fit$Beta[c(5, 17),], tolerance = 1e-3)
      expect_equal(fit_orig$Mu[-c(5, 17),], fit$Mu[-c(5, 17),])
      expect_equal(fit_orig$Mu[c(5, 17),], fit$Mu[c(5, 17),], tolerance = 1e-3)
      expect_equal(fit_orig$overdispersions[-c(5, 17)], fit$overdispersions[-c(5, 17)])
      expect_equal(fit_orig$overdispersions[c(5, 17)], fit$overdispersions[c(5, 17)], tolerance = 0.01)
      expect_equal(fit_orig$overdispersion_shrinkage_list$ql_disp_estimate[-c(5, 17)], fit$overdispersion_shrinkage_list$ql_disp_estimate[-c(5, 17)])
      expect_equal(fit_orig$overdispersion_shrinkage_list$ql_disp_estimate[c(5, 17)], fit$overdispersion_shrinkage_list$ql_disp_estimate[c(5, 17)], tolerance = 0.01)

      expect_equal(fit_orig[-which(names(fit_orig) %in% c("Beta", "overdispersions", "overdispersion_shrinkage_list", "Mu"))],
                   fit[-which(names(fit_orig) %in% c("Beta", "overdispersions", "overdispersion_shrinkage_list", "Mu"))])
  })



  testthat::with_mock(
    # Simulated failure to converge
    fitBeta_fisher_scoring =  function(Y, model_matrix, exp_offset_matrix, thetas, beta_matSEXP, ridge_penalty_nl, tolerance, max_rel_mu_change, max_iter) {
      res <- .Call(`_glmGamPoi_fitBeta_fisher_scoring`, Y, model_matrix, exp_offset_matrix, thetas, beta_matSEXP, ridge_penalty_nl, tolerance, max_rel_mu_change, max_iter)
      if(! any(res$iter == 0)){
        res$beta_mat[c(5, 17), ] <- NA
        res$iter[c(5, 17)] <- 1000
        res$deviance[c(5, 17)] <- NaN
      }
      res
    }, # Make sure that optim doesn't recover the correct result
      estimate_betas_fisher_scoring = {
        ebfs <- estimate_betas_fisher_scoring
        args <- formals(ebfs)
        args[["try_recovering_convergence_problems"]] <- FALSE
        formals(ebfs) <- args
        ebfs
      },
    {
      # The function isn't crashed by internal NA's
      expect_warning(fit <- glm_gp(Y, design = ~ cont + cont2, size_factors = sf))
      expect_equal(sum(is.na(fit$Mu)), 2 * n_cells)
      expect_equal(unname(which(is.na(fit$overdispersions))), c(5, 17))
      expect_equal(unname(which(is.na(fit$overdispersion_shrinkage_list$ql_disp_shrunken))), c(5, 17))
      expect_equal(unname(which(is.na(fit$Beta[, "Intercept"]))), c(5, 17))
      expect_equal(unname(which(is.na(fit$Beta[, "cont"]))), c(5, 17))

      expect_equal(fit_orig$Beta[-c(5, 17),], fit$Beta[-c(5, 17),], tolerance = 1e-4)
      expect_equal(fit_orig$Mu[-c(5, 17),], fit$Mu[-c(5, 17),], tolerance = 1e-4)
      expect_equal(fit_orig$overdispersions[-c(5, 17)], fit$overdispersions[-c(5, 17)])
      expect_equal(fit_orig$overdispersion_shrinkage_list$ql_disp_estimate[-c(5, 17)],
                   fit$overdispersion_shrinkage_list$ql_disp_estimate[-c(5, 17)], tolerance = 1e-3)

      expect_silent(res <- test_de(fit, reduced_design = ~ 1))
      expect_true(is.na(res$pval[5]))
      expect_false(is.na(res$pval[4]))

      expect_silent(res <- test_de(fit, contrast = cont))
      expect_true(is.na(res$pval[5]))
      expect_false(is.na(res$pval[4]))
    })
})



test_that("glm_gp can handle full offset matrices correctly", {


  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))

  sf <- rnorm(n = 50, mean = 1, sd = 0.1)
  mult <- rlnorm(n = 10, meanlog = log(10), sdlog = 1)
  Offset <- log(matrix(mult, ncol = 1) %*% matrix(sf, nrow = 1))

  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = data)
  fit <- glm_gp(se, design = ~ fav_food + city + age, size_factors = sf)
  fit_offset <- glm_gp(se, design = ~ fav_food + city + age, offset = Offset, size_factors = FALSE)

  expect_equal(fit$Beta[,-1], fit_offset$Beta[,-1], tolerance = 1e-3)
  expect_equal(fit$Beta[,1], fit_offset$Beta[,1] + log(mult), tolerance = 1e-3)
})



test_that("glm_gp can handle gigantic counts", {
  skip_on_os("windows", "test fails on windows")
  # Integers with values larger than 2^31-1 = 2,147,483,647 = 2.1 * 10^9 cannot be represented by ints
  mat <- matrix(rpois(n = 5 * 100, lambda = 3 * 10^9), nrow = 5, ncol = 100)
  expect_true(any(mat > .Machine$integer.max))
  rand <- rnorm(100)
  fit <- glmGamPoi::glm_gp(mat, design = ~ rand)
  expect_equal(exp(fit$Beta[,1]), rep(3e9, times = 5), tolerance = 0.1)
  expect_equal(fit$Beta[,2], rep(0, times = 5), tolerance = 1e-4)
})



test_that("providing the col_data explicitly and implicitly doesn't cause issues", {
  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 10 * 50, mu = 3, size = 1/3.1), nrow = 10, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(10))
  colnames(Y) <- paste0("person_", seq_len(50))

  se <- SummarizedExperiment::SummarizedExperiment(Y, colData = data)

  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
  fit_se <- glm_gp(se, design = ~ fav_food + city + age, col_data = SummarizedExperiment::colData(se))
  expect_equal(fit, fit_se)
  expect_equal(colnames(SummarizedExperiment::colData(fit_se$data)), c("fav_food", "city", "age"))

  cd_mod <- SummarizedExperiment::colData(se)
  cd_mod$age[1] <- 42
  expect_warning(
    fit_se2 <- glm_gp(se, design = ~ fav_food + city + age, col_data = cd_mod)
  )
})

