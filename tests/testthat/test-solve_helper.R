set.seed(1)

test_that("solve solve_lm_for_A works", {
  A <- matrix(rnorm(n = 1  * 5), nrow = 1, ncol = 5)
  B <- matrix(rnorm(n = 5 * 20), nrow = 5, ncol = 20)
  y <- A %*% B + rnorm(20, sd = 0.1)

  exact <- solve_lm_for_A(y, B)
  iterative <- optim(par = rep(1, times = 5), function(par){
    sum((y - par %*% B)^2)
  }, method = "BFGS")

  expect_equal(c(exact), iterative$par, tolerance = 1e-6)
  expect_equal(sum((y - exact %*% B)^2), iterative$value)

  A <- matrix(rnorm(n = 3 * 5), nrow = 3, ncol = 5)
  Y <- A %*% B
  expect_equal(A, solve_lm_for_A(Y, B))

  A <- matrix(rnorm(n = 0), nrow = 0, ncol = 5)
  Y <- A %*% B
  expect_equal(A, solve_lm_for_A(Y, B))
})



test_that("solve solve_lm_for_B works", {
  A <- matrix(rnorm(n = 20  * 5), nrow = 20, ncol = 5)
  B <- matrix(rnorm(n = 5 * 1), nrow = 5, ncol = 1)
  y <- A %*% B + rnorm(20, sd = 0.1)

  exact <- solve_lm_for_B(y, A)
  iterative <- optim(par = rep(1, times = 5), function(par){
    sum((y - A %*% par)^2)
  }, method = "BFGS")

  expect_equal(c(exact), iterative$par, tolerance = 1e-6)
  expect_equal(sum((y - A %*% exact)^2), iterative$value)

  B <- matrix(rnorm(n = 5 * 3), nrow = 5, ncol = 3)
  Y <- A %*% B
  expect_equal(B, solve_lm_for_B(Y, A))

  B <- matrix(rnorm(n = 0), nrow = 5, ncol = 0)
  Y <- A %*% B
  expect_equal(B, solve_lm_for_B(Y, A))
})



test_that("weighted solving works", {
  # Check for A
  A <- matrix(rnorm(n = 1  * 5), nrow = 1, ncol = 5)
  B <- matrix(rnorm(n = 5 * 20), nrow = 5, ncol = 20)
  weights <- runif(ncol(B))
  y <- A %*% B + rnorm(20, sd = 0.1)

  exact <- solve_lm_for_A(y, B, w = weights)
  iterative <- optim(par = rep(1, times = 5), function(par){
    sum(weights * (y - par %*% B)^2)
  }, method = "BFGS")

  expect_equal(c(exact), iterative$par, tolerance = 1e-6)
  expect_equal(sum(weights * (y - exact %*% B)^2), iterative$value)


  # Check for B
  A <- matrix(rnorm(n = 20  * 5), nrow = 20, ncol = 5)
  B <- matrix(rnorm(n = 5 * 1), nrow = 5, ncol = 1)
  weights <- runif(nrow(A))
  y <- A %*% B + rnorm(20, sd = 0.1)

  exact <- solve_lm_for_B(y, A, w = weights)
  iterative <- optim(par = rep(1, times = 5), function(par){
    sum(weights * (y - A %*% par)^2)
  }, method = "BFGS")

  expect_equal(c(exact), iterative$par, tolerance = 1e-4)
  expect_equal(sum(weights * (y - A %*% exact)^2), iterative$value)

 })




