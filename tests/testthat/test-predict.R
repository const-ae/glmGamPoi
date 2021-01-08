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



test_that("predict works with vector design", {
  set.seed(1)
  y <- rnbinom(n = 100, mu = 15, size  = 1/0.8)
  group <- sample(LETTERS[1:3], size = 100, replace = TRUE)

  fit_glmGamPoi <- glm_gp(y, design = group)
  pred <- predict(fit_glmGamPoi, newdata = c("A", "A", "B"))
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[c(1, 1, 2)]))
  form <- fit_glmGamPoi$design_formula
  expect_equal(attr(form, "xlevels"), list(x_ = LETTERS[1:3]))

  fit_glmGamPoi <- glm_gp(y, design = group, reference_level = "B")
  pred <- predict(fit_glmGamPoi, newdata = c("A", "A", "B"))
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[1] + c(fit_glmGamPoi$Beta[c(2,2)], 0)))
  form <- fit_glmGamPoi$design_formula
  expect_equal(attr(form, "xlevels"), list(x_ = c("B", "A", "C")))

  int_group <- pmin(rpois(n = 100, lambda = 2), 4) - 2
  fit_glmGamPoi <- glm_gp(y, design = int_group)
  pred <- predict(fit_glmGamPoi, newdata = 0)
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[3]))
  form <- fit_glmGamPoi$design_formula
  expect_equal(attr(form, "xlevels"), list(x_ = as.character(-2:2)))


  fit_glmGamPoi <- glm_gp(y, design = group)
  group2 <- factor(c("C"), levels = c("C", "asdf"))
  pred <- predict(fit_glmGamPoi, newdata = group2[1])
  expect_equal(predict(fit_glmGamPoi, newdata = "C"),  predict(fit_glmGamPoi, newdata = group2[1]))
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[3]))


  group3 <- factor(group, c("C", "A", "B", "D"))
  # No error although "D" is not in data!
  fit_glmGamPoi <- glm_gp(y, design = group3)
  pred <- predict(fit_glmGamPoi, newdata = "C")
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[1]))



  group4 <- rep("A", 100)
  fit_glmGamPoi <- glm_gp(y, design = group4)
  # Maybe this should throw an error because C and G are not in training
  pred <- predict(fit_glmGamPoi, newdata = c("C", "G"))
  expect_equal(drop(pred), unname(fit_glmGamPoi$Beta[c(1,1)]))


})

test_that("predict provides helpful error messages", {
  set.seed(2)
  y <- rnbinom(n = 100, mu = 15, size  = 1/0.8)
  group <- sample(LETTERS[1:3], size = 100, replace = TRUE)

  fit_glmGamPoi <- glm_gp(y, design = group)
  expect_error(predict(fit_glmGamPoi, newdata =  c("A", "B", "F")))


  fit_glmGamPoi <- glm_gp(y, design = ~ group)
  expect_error(predict(fit_glmGamPoi, newdata = data.frame(group = c("A", "B", "F", "G"))))

  cont <- rnorm(100)
  fit_glmGamPoi <- glm_gp(y, design = ~ group + cont)
  expect_error(predict(fit_glmGamPoi, newdata = data.frame(groups = c("A", "B", "F"))))


})



test_that("predict can handle hdf5 arrays", {

  Y <- matrix(rnbinom(n = 100 * 5, mu = 5, size  = 1/0.8), nrow = 5, ncol = 100)
  Y_hdf5 <- HDF5Array::writeHDF5Array(Y)

  col_df <- data.frame(cont = rnorm(100),
                       cont2 = rnorm(100))

  fit_in_memory <- glm_gp(Y, design = ~ cont + cont2, col_data = col_df)
  fit_on_disk <- glm_gp(Y_hdf5, design = ~ cont + cont2, col_data = col_df)

  expect_equal(fit_in_memory[!names(fit_in_memory) %in% c("Mu", "data", "Offset")], fit_on_disk[!names(fit_on_disk) %in% c("Mu", "data", "Offset")])
  expect_equal(c(fit_in_memory$Mu), c(fit_on_disk$Mu))
  expect_s4_class(fit_on_disk$Mu, "DelayedArray")
  expect_equal(c(assay(fit_in_memory$data)), c(assay(fit_on_disk$data)))
  expect_equal(class(fit_in_memory$data), class(fit_on_disk$data))
  expect_s4_class(assay(fit_on_disk$data), "DelayedArray")
  expect_equal(c(fit_in_memory$Offset), c(fit_on_disk$Offset))
  expect_s4_class(fit_on_disk$Offset, "DelayedArray")

  new_data <- data.frame(cont = c(4, 6, 9), cont2 = c(3, 5, 9))
  pred_in_memory <- predict(fit_in_memory, newdata = new_data)
  pred_in_memory2 <- predict(fit_on_disk, newdata = new_data)
  pred_on_disk <- predict(fit_on_disk, newdata = new_data, on_disk = TRUE)

  expect_equal(pred_in_memory, pred_in_memory2)
  expect_s4_class(pred_on_disk, "DelayedArray")
  expect_equal(pred_in_memory, as.matrix(pred_on_disk))



  pred_in_memory <- predict(fit_in_memory, newdata = new_data, se.fit = TRUE)
  pred_in_memory2 <- predict(fit_on_disk, newdata = new_data, se.fit = TRUE)
  pred_on_disk <- predict(fit_on_disk, newdata = new_data, se.fit = TRUE, on_disk = TRUE)

  expect_equal(pred_in_memory, pred_in_memory2)
  expect_s4_class(pred_on_disk$fit, "DelayedArray")
  expect_s4_class(pred_on_disk$se.fit, "DelayedArray")
  expect_equal(pred_in_memory[1:2], lapply(pred_on_disk[1:2], as.matrix))


})



test_that("predict can handle ill-formed weighted design", {

  y <- rep(0, 50)
  cont <- rnorm(n = 50, mean = 100)
  fit <- glm_gp(y ~ cont)
  pred <- predict(fit, se.fit = TRUE, type = "response")
  expect_true(all(pred$fit < 1e-7))
  expect_true(all(pred$se.fit < 1e-5))

  group <- sample(LETTERS[1:5], 50, replace = TRUE)
  fit <- glm_gp(y ~ group)
  pred <- predict(fit, se.fit = TRUE, type = "response")
  expect_true(all(pred$fit == 0))
  expect_true(all(is.na(pred$se.fit)))

})



