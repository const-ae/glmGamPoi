test_that("division with fix for zero by zero division works", {

  expect_equal(div_zbz_int(0, 0), 0)
  expect_equal(div_zbz_dbl(0, 0), 0)
  expect_equal(div_zbz_dbl(pi, exp(1)), pi / exp(1))
  expect_equal(div_zbz_int(3, 4), 3/4)


  # a <- rpois(n = 1e6, lambda = 0.1)
  # b <- rpois(n = 1e6, lambda = 0.1)
  # mean(a == 0)
  # mean(b == 0)
  # bench::mark(div = a / b,
  #             div_zbz = div_zbz_int(a, b),
  #             check = FALSE)
  # # With -g -O2 optimization compilation
  # # expression       min      median
  # #   1 div          2.42ms   2.93ms
  # #   2 div_zbz      3.19ms   4.25ms
  #
  # a <- a * 1.0
  # b <- b * 1.0
  # bench::mark(div = a / b,
  #             div_zbz = div_zbz_dbl(a, b),
  #             check = FALSE)
  # # With -g -O2 optimization compilation
  # # expression       min      median
  # # div              1.23ms   1.66ms
  # # div_zbz          1.46ms   2.21ms
})


test_that("div_zbz_dbl_mat returns a matrix", {
  a <- matrix(rnorm(n = 5 * 3), nrow = 5, ncol = 3)
  b <- matrix(rnorm(n = 5 * 3), nrow = 5, ncol = 3)

  expect_equal(div_zbz_dbl_mat(a, b), a / b)
})
