
# This function was copied from proDA
parse_contrast <- function(contrast, levels, direct_call=TRUE){

  env <- if(direct_call){
    environment()
  }else{
    parent.frame()
  }

  if(missing(contrast)){
    stop("No contrast argument was provided! The option is any linear combination of:\n",
         paste0(levels, collapse = ", "))
  }
  cnt_capture <- substitute(contrast, env = env)


  stopifnot(! is.null(levels))
  if(is.factor(levels)){
    levels <- levels(levels)
  }else if(! is.character(levels)){
    stop("levels must be either a character vector or a factor")
  }

  indicators <- diag(nrow=length(levels))
  rownames(indicators) <- levels
  colnames(indicators) <- levels

  level_environment <- new.env(parent = parent.frame(n = 1 + (!direct_call)))

  for(lvl in levels){
    ind <- indicators[, lvl]
    names(ind) <- levels
    assign(lvl, ind, level_environment)
  }
  tryCatch({
    res <- eval(cnt_capture, envir= level_environment)
    if(! is.numeric(res)){
      if(is.character(res)){
        # If contrast was a string, eval will just spit it out the same way
        res <- eval(parse(text = res), envir= level_environment)
      }
    }
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Object '", match[2], "' not found. Allowed variables in contrast are:\n",
           paste0(levels, collapse = ", "))
    }else{
      stop(e$message)
    }
  })
  res
}




