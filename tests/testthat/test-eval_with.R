


test_that("eval_with works in simple cases", {

  data <- data.frame(x = 1:3, a = LETTERS[1:3], stringsAsFactors = FALSE)

  expect_equal(eval_with(x * 3, data), data$x * 3)
  expect_equal(eval_with(a == "b", data), data$a == "b")
  tmp <- 5
  expect_equal(eval_with(x * tmp, data), data$x * tmp)
  expect_equal(eval_with("x * tmp", data), data$x * tmp)
  string <- "x * tmp"
  expect_equal(eval_with(string, data), data$x * tmp)
  expect_equal(eval_with(a[1:2], data), data$a[1:2])
  expect_error(eval_with(a[1], data))

})



