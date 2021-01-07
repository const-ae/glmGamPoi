
#' @export
predict.glmGamPoi <- function(object, newdata = NULL,
                              type = c("link", "response"),
                              se.fit = FALSE,
                              offset = mean(object$Offset),
                              on_disk = NULL,...){

  type <- match.arg(type, c("link", "response"))
  if(is.null(newdata)){
    # Easy, just return mu
    Mu <- object$Mu
    design_matrix <- object$model_matrix
  }else{
    # Do something with newdata
    if(is.matrix(newdata)){
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

    Mu <- calculate_mu(object$Beta, design_matrix, offset_matrix)
    rownames(Mu) <- rownames(object$Beta)
    colnames(Mu) <- mu_colnames
  }

  if(type == "response"){
    fit <- Mu
  }else if(type == "link"){
    fit <- log(Mu)
  }


  if(se.fit){
    p_idxs <- seq_len(ncol(object$model_matrix))
    # This could be adapted to the quasi-GamPoi value
    scale <- 1
    se_fit <- t(vapply(seq_len(nrow(fit)), function(gene_idx){
      disp <- object$overdispersions[gene_idx]
      mu <- object$Mu[gene_idx, ]
      weights <- mu / (1 + mu * disp)
      weighted_Design <-  object$model_matrix * sqrt(weights)
      R <- qr.R(qr(weighted_Design))[p_idxs, p_idxs, drop=FALSE]
      XRinv <- design_matrix %*% qr.solve(R)
      sqrt(matrixStats::rowSums2(XRinv^2))
    }, FUN.VALUE = rep(0.0, ncol(fit))))

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
    stop("Cannot handle offset of class '", class(data), "'.",
         "It must be either a scalar/numeric vector or dense matrix ",
         "object (i.e., a base matrix or DelayedArray).")
  }
  offset_mat
}


