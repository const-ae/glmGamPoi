
eval_with <- function(statement, data){
  stmnt_capture <- substitute(statement)
  eval_with_q(stmnt_capture, data = data, env = parent.frame())
}

eval_with_q <- function(statement, data, env = parent.frame()){
  statement_envir <- new.env(parent = env)
  tryCatch({
    if(is.character(statement)){
      res <- eval(parse(text = statement), envir= as.list(data), enclos = statement_envir)
    }else{
      res <- eval(statement, envir = as.list(data), enclos = statement_envir)
    }
  }, error = function(e){
    # Try to extract text from error message
    match <- regmatches(e$message, regexec("object '(.+)' not found", e$message))[[1]]
    if(length(match) == 2){
      stop("Object '", match[2], "' not found. Allowed variables are:\n",
           paste0(names(data), collapse = ", "))
    }else{
      stop(e$message)
    }
  })
  res
}

