test_that("glm_gp works for simple cases", {
  y <- rpois(n = 10, lambda = 3)
  expect_warning({
    tmp <- glm_gp(y, design = ~ 1)
  })
  tmp <- glm_gp(y, design = ~ 1, size_factors = FALSE)
  expect_equal(tmp$size_factors, rep(1, times = 10))
  expect_equal(c(tmp$model_matrix), rep(1, times = 10))
  expect_equal(c(tmp$Mu_est), rep(mean(y), times = 10))
  expect_equal(c(tmp$Beta_est), log(mean(y)))
})


test_that("glm_gp can handle complicated formulas", {

  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                     age = rnorm(n = 50, mean = 40, sd = 15))
  Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
  rownames(Y) <- paste0("gene_", seq_len(100))
  colnames(Y) <- paste0("person_", seq_len(50))

  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
  expect_equal(colnames(fit$Beta_est), colnames(fit$model_matrix))
  expect_equal(rownames(fit$Beta_est), rownames(Y))
  expect_equal(colnames(fit$Mu_est), colnames(Y))
  expect_equal(rownames(fit$Mu_est), rownames(Y))
  expect_equal(rownames(fit$model_matrix), colnames(Y))
  expect_equal(names(fit$overdispersions), rownames(Y))
  expect_equal(names(fit$size_factors), colnames(Y))
})


test_that("glm_gp can handle no-row input", {

  Y <- matrix(numeric(0), nrow = 0, ncol = 10)
  tmp <- glm_gp(Y, size_factors = FALSE)
  expect_equal(dim(tmp$Beta_est), c(0, 1))
  expect_equal(dim(tmp$Mu_est), c(0, 10))
  expect_equal(tmp$size_factors, rep(1, times = 10))
  expect_equal(c(tmp$model_matrix), rep(1, times = 10))


})

test_that("glm_gp can handle no-col input", {

  Y <- matrix(numeric(0), nrow = 3, ncol = 0)
  tmp <- glm_gp(Y, size_factors = FALSE)
  expect_equal(dim(tmp$Beta_est), c(3, 1))
  expect_equal(c(tmp$Beta_est), rep(-Inf, times = 3))
  expect_equal(dim(tmp$Mu_est), c(3, 0))
  expect_equal(tmp$size_factors, rep(1, times = 0))
  expect_equal(dim(tmp$model_matrix), c(0, 1))
  expect_equal(c(tmp$model_matrix), rep(1, times = 0))


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
