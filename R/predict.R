
#' Predict 'link' or 'response' values for Gamma-Poisson GLMs
#'
#' Predict `mu` (i.e., `type = "response"`) or `log(mu)` (i.e., `type = "link"`)
#' from a 'glmGamPoi' fit (created by `glm_gp(...)`) with the corresponding
#' estimate of the standard error. If `newdata` is `NULL`, `mu` is returned
#' for the original input data.
#'
#' @param object a `glmGamPoi` fit object (produced by `glm_gp()`).
#' @param newdata a specification of the new data for which the
#'   expression for each gene is predicted. `newdata` should be a
#'   \describe{
#'     \item{data.frame}{if the original fit was specified with
#'       a `formula`, provide a `data.frame` with one column for
#'       each variable in the formula.
#'       For example, if `glm_gp(se, design = ~ age + batch + treatment)`,
#'       then the data.frame needs a `age`, `batch`, and `treatment` column that
#'       contain the same data types as the original fit.}
#'    \item{vector}{if the original fit was specified using a vector, you
#'       need to again provide a vector with the same format.}
#'    \item{matrix}{if `newdata` is a matrix, it is applied directly
#'      as `Mu <- exp(object$Beta %*% t(newdata) + object$offset_matrix)`.
#'      So make sure, that it is constructed correctly.}
#'    \item{NULL}{if `newdata` is `NULL`, the predicted values for the original
#'      input data are returned.}
#'   }
#' @param type either 'link' or 'response'. The default is 'link', which returns the predicted values
#'   before the link function (`exp()`) is applied. Thus, the values can be
#'   positive and negative numbers. However, often the predicted values are easier to interpret
#'   **after** the link function is applied (i.e., type = "response"), because then the
#'   values are on the same scale as the original counts.
#' @param se.fit boolean that indicates if in addition to the mean the standard error of the
#'   mean is returned.
#' @param offset count models (in particular for sequencing experiments) usually have a sample
#'   specific size factor (`offset = log(size factor)`). It defines how big we expect the predicted
#'   results are. If `newdata` is `NULL`, the `offset` is ignored, because the `predict()` returns
#'   a result based on the pre-calculated `object$Mu`. If `newdata` is not `NULL`, by default the
#'   `offset` is `mean(object$Offset)`, which puts the in the same size as the average sample.
#' @param on_disk a boolean that indicates if the results are `HDF5Matrix`'s from the `HDF5Array`
#'   package. If `newdata` is `NULL`, `on_disk` is ignored. Otherwise, if `on_disk = NULL`,
#'   the result is calculated on disk depending if `offset` is stored on disk.
#' @param verbose a boolean that indicates if information about the individual steps are
#'   printed while predicting. Default: `FALSE`.
#' @param ... currently ignored.
#'
#' @details
#'   For `se.fit = TRUE`, the function sticks very close to the behavior of `stats::predict.glm()` for
#'   fits from `MASS::glm.nb()`.
#'
#' @return
#'   If `se.fit == FALSE`, a matrix with dimensions `nrow(object$data) x nrow(newdata)`. \cr
#'   If `se.fit == TRUE`, a list with three entries
#'   \describe{
#'     \item{fit}{the predicted values as a matrix with dimensions `nrow(object$data) x nrow(newdata)`.
#'       This is what would be returned if `se.fit == FALSE`.}
#'     \item{se.fit}{the associated standard errors for each `fit`. Also a matrix with
#'        dimensions `nrow(object$data) x nrow(newdata)`.}
#'     \item{residual.scale}{Currently fixed to 1. In the future, this might become the values from
#'        `object$overdispersion_shrinkage_list$ql_disp_shrunken`.}
#'   }
#'
#' @seealso [stats::predict.lm()] and [stats::predict.glm()]
#'
#' @examples
#'
#'  set.seed(1)
#'  # The simplest example
#'  y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
#'  fit <- glm_gp(y, size_factors = FALSE)
#'  predict(fit, type = "response")
#'  predict(fit, type = "link", se.fit = TRUE)
#'
#'
#'  # Fitting a whole matrix
#'  model_matrix <- cbind(1, rnorm(5))
#'  true_Beta <- cbind(rnorm(n = 30), rnorm(n = 30, mean = 3))
#'  sf <- exp(rnorm(n = 5, mean = 0.7))
#'  model_matrix
#'  Y <- matrix(rnbinom(n = 30 * 5, mu = sf * exp(true_Beta %*% t(model_matrix)), size = 1/2.4),
#'              nrow = 30, ncol = 5)
#'  fit <- glm_gp(Y, design = model_matrix, size_factors = sf, verbose = TRUE)
#'
#'  head(predict(fit, type = "response"))
#'  pred <- predict(fit, type = "link", se.fit = TRUE, verbose = TRUE)
#'  head(pred$fit)
#'  head(pred$se.fit)
#'
#'
#'  # Fitting a model with covariates
#'  data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
#'                     city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
#'                     age = rnorm(n = 50, mean = 40, sd = 15))
#'  Y <- matrix(rnbinom(n = 4 * 50, mu = 3, size = 1/3.1), nrow = 4, ncol = 50)
#'  fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
#'  predict(fit)[, 1:3]
#'
#'  nd <- data.frame(fav_food = "banana", city = "paris", age = 29)
#'  predict(fit, newdata = nd)
#'
#'  nd <- data.frame(fav_food = "banana", city = "paris", age = 29:40)
#'  predict(fit, newdata = nd, se.fit = TRUE, type = "response")
#'
#'
#' @export
predict.glmGamPoi <- function(object, newdata = NULL,
                              type = c("link", "response"),
                              se.fit = FALSE,
                              offset = mean(object$Offset),
                              on_disk = NULL, verbose = FALSE,
                              ...){

  type <- match.arg(type, c("link", "response"))
  if(is.null(newdata)){
    # Easy, just return mu
    if(verbose) message("'newdata' is NULL, use 'Mu = object$Mu'")
    Mu <- object$Mu
    design_matrix <- object$model_matrix
  }else{
    # Do something with newdata
    if(is.matrix(newdata)){
      if(verbose) message("'newdata' is a matrix, set 'design_matrix = newdata'")
      design_matrix <- newdata
      mu_colnames <- rownames(newdata)
    }else  if((is.vector(newdata) || is.factor(newdata))){
      cf <- attr(object$design_formula, "constructed_from")
      if(is.null(cf) || cf != "vector"){
        stop("Original model was constructed from a ", cf, ". Please, ",
             "provide 'newdata' as a dataframe and not as a vector.",
             call. = FALSE)
      }
      if(is.numeric(newdata)){
        newdata <- as.character(newdata)
      }
      if(verbose) message("'newdata' is a vector, construct 'design_matrix' from 'newdata' and 'object$design_formula'")
      mu_colnames <- names(newdata)
      newdata <- data.frame(x_ = newdata, stringsAsFactors = FALSE)
      design_matrix <- make_model_matrix_for_predict(object, newdata)
    }else if(is.data.frame(newdata)){
      cf <- attr(object$design_formula, "constructed_from")
      if(is.null(cf) || cf != "formula"){
        stop("The original model was constructed from a ", cf, ". Please, ",
             "provide 'newdata' as a vector and not as a data.frame.",
             call. = FALSE)
      }
      if(verbose) message("'newdata' is a data.frame, construct 'design_matrix' from 'newdata' and 'object$design_formula'")
      design_matrix <- make_model_matrix_for_predict(object, newdata)
      mu_colnames <- rownames(newdata)
    }else{
      stop("Don't know how to handle newdata of class ",
           paste0(class(newdata), collapse = ", "), ". Please provide ",
           "a 'data.frame' if the original model was constructed from a 'formula' ",
           "or a 'vector' if the original model was constructed from a 'vector'.",
           call. = FALSE)
    }

    offset_matrix <- handle_offset_param_for_predict(offset, nrow = nrow(object$Beta),
                                    ncol = nrow(design_matrix), on_disk = on_disk)
    if(verbose) message("Calculate 'Mu = exp(object$Beta %*% t(design_matrix) + Offset)'")
    Mu <- calculate_mu(object$Beta, design_matrix, offset_matrix)
    rownames(Mu) <- rownames(object$Beta)
    colnames(Mu) <- mu_colnames
  }

  if(type == "response"){
    if(verbose) message("'type = \"response\"', return Mu")
    fit <- Mu
  }else if(type == "link"){
    if(verbose) message("'type = \"link\"', return 'fit = log(Mu)'")
    fit <- log(Mu)
  }


  if(se.fit){
    if(verbose) message("'se.fit' is TRUE, calculate the standard error for each Mu estimate")
    p_idxs <- seq_len(ncol(object$model_matrix))
    # This could (should?) be adapted to the quasi-GamPoi value
    scale <- 1
    wrtm <- if( is(fit, "DelayedMatrix") && is(DelayedArray::seed(fit), "HDF5ArraySeed")){
      write_rows_to_hdf5_matrix
    }else{
      write_rows_to_matrix
    }
    # This could be made more efficient by batching the reads
    # of object$Mu
    se_fit <- wrtm(nrow(fit), ncol(fit), object$Mu, function(gene_idx, mu){
      disp <- object$overdispersions[gene_idx]
      weights <- mu / (1 + mu * disp)
      weighted_Design <-  object$model_matrix * sqrt(weights)
      tryCatch({
        R <- qr.R(qr(weighted_Design))[p_idxs, p_idxs, drop=FALSE]
        XRinv <- design_matrix %*% qr.solve(R)
        sqrt(matrixStats::rowSums2(XRinv^2))
      }, error = function(err){
        # For example R is singular
        # This can happen if all mu == 0
        rep(NA_real_, nrow(design_matrix))
      })
    })

    if(type == "response"){
      se_fit <- se_fit * Mu
    }
    dimnames(se_fit) <- dimnames(fit)
    list(fit = fit, se.fit = se_fit, residual.scale = scale)
  }else{
    fit
  }
}



make_model_matrix_for_predict <- function(object, newdata){
  stopifnot("glmGamPoi" %in% class(object))
  form <- object$design_formula
  if(is.null(form)){
    stop("predict.glmGamPoi was called with 'newdata' that is not a matrix. ",
         "However, the original model was constructed directly from a matrix, so there is ",
         "no formula to convert the newdata to an appropriate design matrix.",
         call. = FALSE)
  }
  if(! inherits(form, "terms") ||
      is.null(attr(form, "constructed_from", exact = TRUE))){
    stop("object$design_formula must be a terms object and have ",
         " a 'constructed_from' attribute.")
  }

  form <- delete.response(form)
  tryCatch({
    mf <- model.frame(form, newdata, xlev = attr(form, "xlevels"))
  }, error = function(e){
    match_new_factor <- regmatches(e$message, regexec("factor (.+) has new levels (.+)", e$message))[[1]]
    match_obj_not_found <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match_new_factor) == 3){
      if(match_new_factor[2] == "x_"){
        stop("while constructing design matrix from 'newdata'.\n",
             "'newdata' contained values (", match_new_factor[3], ") which were not in the training data.\n",
             "The training data contained values:\n",
             paste0(attr(form, "xlevels")[["x_"]], collapse = ", "), call. = FALSE)
      }else{
        stop("while constructing design matrix from 'newdata'\n",
             "Column '", match_new_factor[2], "' of the 'newdata' contained values (", match_new_factor[3], ") which were not in the training data.\n",
             "The training data contained values:\n",
             paste0(attr(form, "xlevels")[[match_new_factor[2]]], collapse = ", "), call. = FALSE)
      }
    }else if(length(match_obj_not_found) == 2){
      stop("while constructing design matrix from 'newdata'\n",
           "The column '", match_obj_not_found[2],"' was absent from 'newdata'", call. = FALSE)
    }else{
      stop(e$message)
    }
  })
  mm <- model.matrix(form, mf, constrasts.arg = attr(object$model_matrix, "contrasts", exact = TRUE))

  mm
}



handle_offset_param_for_predict <- function(offset, nrow, ncol, on_disk){
  if(is.numeric(offset)){
    stopifnot(length(offset) == 1 || length(offset) == ncol)
    offset <- matrix(offset, nrow = nrow, ncol = ncol, byrow = TRUE)
  }
  # Check that offset is correctly sized
  stopifnot(nrow(offset) == nrow && ncol(offset) == ncol)

  if(is.matrix(offset)){
    if(is.null(on_disk) || isFALSE(on_disk)){
      offset_mat <- offset
    }else if(isTRUE(on_disk)){
      offset_mat <- HDF5Array::writeHDF5Array(offset)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else if(is(offset, "DelayedArray")){
    if(is.null(on_disk) || isTRUE(on_disk)){
      offset_mat <- offset
    }else if(isFALSE(on_disk)){
      warning("offset was provided as a DelayedArray, but 'on_disk = FALSE'. ",
              "I will thus realize the full matrix in memory.")
      offset_mat <- as.matrix(offset)
    }else{
      stop("Illegal argument type for on_disk. Can only handle 'NULL', 'TRUE', or 'FALSE'")
    }
  }else{
    stop("Cannot handle offset of class '", class(offset), "'.",
         "It must be either a scalar/numeric vector or dense matrix ",
         "object (i.e., a base matrix or DelayedArray).")
  }
  offset_mat
}



write_rows_to_matrix <- function(nrow, ncol, Mu, FUN, ...){
  res <- matrix(nrow = nrow, ncol = ncol)
  for(idx in seq_len(nrow)){
    res[idx, ] <- FUN(idx, Mu[idx, ], ...)
  }
  res
}

write_rows_to_hdf5_matrix <- function(nrow, ncol, Mu, FUN, ...){
  res_sink <- HDF5Array::HDF5RealizationSink(c(nrow, ncol))
  on.exit({
    DelayedArray::close(res_sink)
  }, add = TRUE)
  stopifnot(nrow(Mu) == nrow)
  res_grid <- DelayedArray::rowAutoGrid(res_sink)
  mu_grid <- DelayedArray::rowAutoGrid(Mu)

  res_block_counter <- 1L
  mu_block_counter <- 1L
  res_block <- res_grid[[res_block_counter]]
  mu_block <- mu_grid[[mu_block_counter]]
  res_block_end <- BiocGenerics::end(res_block)[1]
  mu_block_end <- BiocGenerics::end(mu_block)[1]
  Res_subset <- matrix(NA, nrow = BiocGenerics::width(res_block)[1], ncol = ncol)
  Mu_subset <- DelayedArray::read_block(Mu, mu_block)
  Res_subset_idx <- 1L
  Mu_subset_idx <- 1L
  for(gene_idx in seq_len(nrow)){
    if(gene_idx > res_block_end){
      # Increment res_block
      Res_subset_idx <- 1L
      res_block_counter <- res_block_counter + 1L
      res_block <- res_grid[[res_block_counter]]
      res_block_end <- BiocGenerics::end(res_block)[1]
      Res_subset <- matrix(NA, nrow = BiocGenerics::width(res_block)[1], ncol = ncol)
    }
    if(gene_idx > mu_block_end){
      # Increment res_block
      Mu_subset_idx <- 1L
      mu_block_counter <- mu_block_counter + 1L
      mu_block <- mu_grid[[mu_block_counter]]
      mu_block_end <- BiocGenerics::end(mu_block)[1]
      Mu_subset <- DelayedArray::read_block(Mu, mu_block)
    }

    # Call FUN
    Res_subset[Res_subset_idx, ] <- FUN(gene_idx, Mu_subset[Mu_subset_idx, ], ...)

    Res_subset_idx <- Res_subset_idx + 1L
    Mu_subset_idx <- Mu_subset_idx + 1L
    if(gene_idx == res_block_end){
      DelayedArray::write_block(res_sink, res_block, Res_subset)
    }
  }
  as(res_sink, "DelayedArray")
}



