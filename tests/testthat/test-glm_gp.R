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

  expect_equal(fit_in_memory[!names(fit_in_memory) %in% c("Mu", "data")], fit_on_disk[!names(fit_on_disk) %in% c("Mu", "data")])
  expect_equal(c(fit_in_memory$Mu), c(fit_on_disk$Mu))
  expect_s4_class(fit_on_disk$Mu, "DelayedArray")
  expect_equal(c(assay(fit_in_memory$data)), c(assay(fit_on_disk$data)))
  expect_equal(class(fit_in_memory$data), class(fit_on_disk$data))
  expect_s4_class(assay(fit_on_disk$data), "DelayedArray")


  fit_on_disk2 <- glm_gp(Y_hdf5, design = ~ fav_food + city + age, col_data = data)
  expect_equal(fit_on_disk[!names(fit_in_memory) %in% c("Mu", "data")], fit_on_disk2[!names(fit_on_disk2) %in% c("Mu", "data")])
  expect_s4_class(fit_on_disk2$Mu, "DelayedArray")
  expect_s4_class(assay(fit_on_disk2$data), "DelayedArray")

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
  expect_equal(fit_vec, fit_mat)


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

test_that("glm_gp warns about mismatching col_data rownames ", {

  coldata <- data.frame(condition = c(rep("A", 4), rep("B", 3)), stringsAsFactors = FALSE)
  rownames(coldata) <- paste0("sample_", sample(1:7))
  Y <- matrix(numeric(0), ncol = 7, nrow = 0)
  colnames(Y) <- paste0("sample_", 1:7)
  fit <- glm_gp(Y, design = ~ condition, col_data = coldata, size_factors = FALSE)
  expect_equal(c(fit$model_matrix), c(model.matrix(~ condition, coldata)[colnames(Y), ]))
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

