test_that("loc_median_fit works", {
  x <- sample(1:10)
  y <- rnorm(n = 10)
  expect_equal(length(x), length(loc_median_fit(x, y, npoints = 1)))
  expect_equal(length(x), length(loc_median_fit(x, y, npoints = 2)))
  expect_equal(length(x), length(loc_median_fit(x, y, npoints = 10)))
  expect_equal(length(x), length(loc_median_fit(x, y, npoints = 20)))

  r1 <- loc_median_fit(x, y)

  r2 <- loc_median_fit(x, y, fraction = 1, weighted = FALSE)
  expect_equal(r2, rep(median(y), 10))
  r25 <- loc_median_fit(x[-1], y[-1], fraction = 1, weighted = FALSE)
  expect_equal(r25, rep(median(y[-1]), 9))

  r4 <- loc_median_fit(x[1], y[1])
  expect_equal(r4, y[1])

  r5 <- loc_median_fit(numeric(0), numeric(0))
  expect_equal(r5, numeric(0))

  y[c(1, 4)] <- 0

  r6 <- loc_median_fit(x, y, npoints = 1, ignore_zeros = TRUE)
  expect_equal(y == 0, is.na(r6))


  r7 <- loc_median_fit(x, y, npoints = 4, ignore_zeros = TRUE)
  expect_true(! any(is.na(r7)))

})


test_that("loc_median_fit works with subsampling", {
  x <- sample(1:10)
  y <- rnorm(n = 10)
  expect_equal(loc_median_fit(x, y, npoints = 1, sample_fraction = 0.9999),
               loc_median_fit(x, y, npoints = 1, sample_fraction = 1))
  expect_equal(loc_median_fit(x, y, npoints = 2, sample_fraction = 0.9999),
               loc_median_fit(x, y, npoints = 2, sample_fraction = 1))
  expect_equal(loc_median_fit(x, y, npoints = 10, sample_fraction = 0.9999, weighted = TRUE),
               loc_median_fit(x, y, npoints = 10, sample_fraction = 1, weighted = TRUE))
  expect_equal(loc_median_fit(x, y, npoints = 20, sample_fraction = 0.9999),
               loc_median_fit(x, y, npoints = 20, sample_fraction = 1))


  expect_equal(length(loc_median_fit(x, y, fraction = 1, sample_fraction = 1, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, fraction = 1, sample_fraction = 0.9999, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, fraction = 1, sample_fraction = 0.5, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, fraction = 1, sample_fraction = 0.00001, weighted = FALSE)), 10)

  expect_equal(length(loc_median_fit(x, y, npoints = 4, sample_fraction = 1, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, npoints = 4, sample_fraction = 0.9999, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, npoints = 4, sample_fraction = 0.5, weighted = FALSE)), 10)
  expect_equal(length(loc_median_fit(x, y, npoints = 4, sample_fraction = 0.00001, weighted = FALSE)), 10)

  expect_equal(length(loc_median_fit(x[0], y[0], npoints = 4, sample_fraction = 1, weighted = FALSE)), 0)
  expect_equal(length(loc_median_fit(x[0], y[0], npoints = 4, sample_fraction = 0.9999, weighted = FALSE)), 0)
  expect_equal(length(loc_median_fit(x[0], y[0], npoints = 4, sample_fraction = 0.5, weighted = FALSE)), 0)
  expect_equal(length(loc_median_fit(x[0], y[0], npoints = 4, sample_fraction = 0.00001, weighted = FALSE)), 0)
})




