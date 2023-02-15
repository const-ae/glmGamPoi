
parse_contrast <- function(contrast, coefficient_names, formula = NULL) {
  if(missing(contrast)){
    stop("No contrast argument was provided! The option is any linear combination of:\n",
         paste0(coefficient_names, collapse = ", "))
  }
  cnt_capture <- rlang::enquo(contrast)

  stopifnot(! is.null(coefficient_names))
  if(is.factor(coefficient_names)){
    coefficient_names <- levels(coefficient_names)
  }else if(! is.character(coefficient_names)){
    stop("levels must be either a character vector or a factor")
  }

  indicators <- diag(nrow=length(coefficient_names))
  rownames(indicators) <- coefficient_names
  colnames(indicators) <- coefficient_names

  covar_indicators <- list()
  for(lvl in coefficient_names){
    ind <- indicators[, lvl]
    names(ind) <- coefficient_names
    covar_indicators[[lvl]] <- ind
  }
  top <- rlang::new_environment(list(
    cond = function(...){
      .cond(formula, list(...))
    }))
  bottom <- rlang::new_environment(covar_indicators, parent = top)
  data_mask <- rlang::new_data_mask(bottom = bottom, top = top)
  data_mask$.cntrst <- rlang::as_data_pronoun(bottom)
  tryCatch({
    res <- rlang::eval_tidy(cnt_capture, data = data_mask)
    if(! is.numeric(res)){
      if(is.character(res)){
        # If contrast was a string, eval will just spit it out the same way
        res <- rlang::eval_tidy(rlang::parse_expr(res), data = data_mask)
      }
    }
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Object '", match[2], "' not found. Allowed variables in contrast are:\n",
           paste0(coefficient_names, collapse = ", "), call. = FALSE)
    }else{
      stop(e$message)
    }
  })
  res
}



.cond <- function(formula, level_sets = list()){
  if(is.null(formula)){
    stop("You called 'cond()' inside the contrast, however the original model ",
         "was not specified with a formula. Thus 'cond()' doesn't work and you ",
         "need to specify the contrast using the column names of the design matrix.")
  }
  if(is.null(attr(formula, "xlevels"))){
    warning("The formula has no 'xlevels' attribute. This is supicious and might indicate a bug.")
  }
  if(any(names(level_sets) == "")){
    stop("All arguments to 'cond()' must be named.")
  }
  if(any(duplicated(names(level_sets)))){
    stop("All arguments to 'cond()' must be unique.")
  }

  covar <- all.vars(formula)
  new_dat <- as.list(rep(0, length(covar)))
  names(new_dat) <- covar
  xlevels <- attr(formula, "xlevels")
  for(n in names(xlevels)){
    new_dat[[n]] <- factor(xlevels[[n]][1], levels = xlevels[[n]])
  }
  for(n in names(level_sets)){
    if(! n %in% names(new_dat)){
      stop("Setting the level of '", n, "' failed. You can only set the level of the following variables: ", paste0(covar, collapse = ", "))
    }
    if(length(level_sets[[n]]) != 1){
      stop("Each argument to 'cond()' must be length one. '", n, "' has length ", length(level_sets[[n]]))
    }
    if(n %in% names(xlevels)){
      if(! level_sets[[n]] %in% xlevels[[n]]){
        stop("You are trying to set '", n, "=", level_sets[[n]], "'. However only the following values for ", n,
             " are valid: ", paste0(xlevels[[n]], collapse = ", "))
      }
      new_dat[[n]] <- factor(level_sets[[n]], levels = xlevels[[n]])
    }else{
      new_dat[[n]] <- level_sets[[n]]
    }
  }
  res <- drop(model.matrix(formula, new_dat, contrasts.arg = attr(formula, "contrasts")))
  attr(res, "assign") <- NULL
  attr(res, "contrasts") <- NULL
  res
}

