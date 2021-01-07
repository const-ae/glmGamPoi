test_that("predict works for simple cases", {
  set.seed(1)
  y <- rnbinom(n = 100, mu = 15, size  = 1/0.8)
  design <- cbind(1, matrix(rnorm(n = 100 * 4), nrow = 100, ncol = 4))

  # use glm.nb instead of glm(..., family = negative.binomial(theta = 3))
  # because otherwise result isn't tagged with class negbin and thus
  # the dispersion would have to be explcitly set to 1 everywhere
  fit_glm <- MASS::glm.nb(y ~ design - 1)
  fit_glmGamPoi <- glm_gp(y ~ design - 1, overdispersion = 1/fit_glm$theta)

  expect_equal(fit_glm$coefficients, drop(fit_glmGamPoi$Beta), tolerance = 1e-5)
  expect_lte(sum(residuals(fit_glmGamPoi, type = "deviance")^2), sum(residuals(fit_glm, type = "deviance")^2))


  # Compare predict()
  expect_equal(lapply(predict(fit_glm, se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, type = "link", se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, type = "link", se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, type = "response", se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, type = "response", se.fit = TRUE), drop),
               tolerance = 1e-5)

})




test_that("predict works for new data", {
  set.seed(1)
  y <- rnbinom(n = 100, mu = 15, size  = 1/0.8)
  df <- data.frame(group = sample(LETTERS[1:3], size = 100, replace = TRUE),
                   cont = rnorm(100))


  # use glm.nb instead of glm(..., family = negative.binomial(theta = 3))
  # because otherwise result isn't tagged with class negbin and thus
  # the dispersion would have to be explcitly set to 1 everywhere
  fit_glm <- MASS::glm.nb(y ~ group + cont, data = df)
  fit_glmGamPoi <- glm_gp(y, ~ group + cont, col_data = df, overdispersion = 1/fit_glm$theta)

  expect_equal(unname(fit_glm$coefficients), unname(drop(fit_glmGamPoi$Beta)), tolerance = 1e-5)
  # expect_lte(sum(residuals(fit_glmGamPoi, type = "deviance")^2), sum(residuals(fit_glm, type = "deviance")^2), tolerance = 1e-7)

  new_data <- data.frame(group = "B", cont = 3)

  # Compare predict()
  # The unname stuff is necessary, because predict.glm is inconsistent with its results...
  expect_equal(lapply(predict(fit_glm, newdata = new_data, se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, newdata = new_data, se.fit = TRUE), function(t)unname(drop(t))),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, newdata = new_data, type = "link", se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "link", se.fit = TRUE), function(t) unname(drop(t))),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, newdata = new_data, type = "response", se.fit = TRUE), unname),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "response", se.fit = TRUE), function(t)unname(drop(t))),
               tolerance = 1e-5)

  new_data <- df[1:10,,drop=FALSE]
  expect_equal(lapply(predict(fit_glm, newdata = new_data, se.fit = TRUE), identity),
               lapply(predict(fit_glmGamPoi, newdata = new_data, se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, newdata = new_data, type = "link", se.fit = TRUE), identity),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "link", se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(lapply(predict(fit_glm, newdata = new_data, type = "response", se.fit = TRUE), identity),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "response", se.fit = TRUE), drop),
               tolerance = 1e-5)

})




