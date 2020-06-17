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

test_that("parse contrast can be called through multiple functions", {


  c1 <- parse_contrast(A - B, levels = LETTERS[1:2])
  fnc1 <- function(contrast){
      parse_contrast(contrast, levels = LETTERS[1:2], direct_call = FALSE)
  }
  fnc2 <- function(crt){
    cnt_capture <- substitute(crt)
    fnc1(cnt_capture)
  }
  c2 <- fnc2(A - B)
  expect_equal(c1, c2)
  c3 <- fnc2("A - B")
  expect_equal(c1, c3)
  string <- "A - B"
  c4 <- fnc2(string)
  expect_equal(c1, c4)
  c5 <- fnc2(paste0("A","-", "B"))
  expect_equal(c1, c5)

})

