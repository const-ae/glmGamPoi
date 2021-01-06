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


test_that("handle_formula works", {

  df <- data.frame(letters = sample(LETTERS[1:3], size = 100, replace = TRUE),
                   factor = as.factor(sample(letters[1:3], size = 100, replace = TRUE)),
                   cont = rnorm(100))

  handle_design_parameter(~ letters + factor + cont, data = matrix(nrow = 5, ncol = 100), col_data = df, reference_level = NULL)

  convert_chr_vec_to_model_matrix(df$letters, NULL)

  # debugonce(model.frame.default)
  # form <- ~ letters + factor + I(cont * 5)
  # form <- ~ letters*factor
  # model.frame(form, data = df)
  # formula.tools::get.vars(form)
  # te <- terms.formula(form, data = df)
  # attr(te, "term.labels")[attr(te, "order") == 1]
  # head(model.frame(form, df))
  #
  # str(te)
  # attr(te, ".Evironment") <- c()
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
  expect_equal(predict(fit_glm, newdata = new_data, se.fit = TRUE),
               lapply(predict(fit_glmGamPoi, newdata = new_data, se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(predict(fit_glm, newdata = new_data, type = "link", se.fit = TRUE),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "link", se.fit = TRUE), drop),
               tolerance = 1e-5)
  expect_equal(predict(fit_glm, newdata = new_data, type = "response", se.fit = TRUE),
               lapply(predict(fit_glmGamPoi, newdata = new_data, type = "response", se.fit = TRUE), drop),
               tolerance = 1e-5)

})




df <- data.frame(group = sample(LETTERS[1:3], size = 100, replace = TRUE),
                 group2 = sample(LETTERS[1:3], size = 100, replace = TRUE),
                 cont = rnorm(100))
mf <- model.frame(~ group + group2 + cont + group:group2, dat = df)
attr(mf, "terms")
str(mf)
