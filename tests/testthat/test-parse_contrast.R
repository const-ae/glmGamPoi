test_that("parse_contrast works", {

  expect_equal(parse_contrast(A, levels = LETTERS[1:5]),
               c(A=1, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast(A - B, levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast("A - B", levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast(A * 3, levels = LETTERS[1:5]),
               c(A=3, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast("Intercept", levels = c("Intercept")),
               c(Intercept = 1))

})

test_that("Parser contrast can handle reference to object in environment", {
  c1 <- parse_contrast(A - B, levels = LETTERS[1:2])
  c2 <- parse_contrast("A - B", levels = LETTERS[1:2])
  string <- "A - B"
  c3 <- parse_contrast(string, levels = LETTERS[1:2])
  c4 <- parse_contrast(paste0("A", " - ", "B"), levels = LETTERS[1:2])
  expect_equal(c1, c2)
  expect_equal(c1, c3)
  expect_equal(c1, c4)
})


test_that("parse contrast can handle being inside a function", {

  fnc <- function(){
    string <- "A - B"
    parse_contrast(string, levels = LETTERS[1:2])
  }

  c1 <- parse_contrast(A - B, levels = LETTERS[1:2])
  c2 <- fnc()
  expect_equal(c1, c2)
})



test_that("parse contrast throws appropriate errors", {

  expect_error(parse_contrast(A - B, levels = LETTERS[3:4]))
  expect_error(parse_contrast(A - B, levels = 3:4))
  expect_error(parse_contrast(levels = LETTERS[1:2]))

})



