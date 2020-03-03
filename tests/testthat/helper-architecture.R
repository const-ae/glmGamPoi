

is_macos <- function(){
  Sys.info()["sysname"] == "Darwin"
}

is_windows <- function(){
  Sys.info()["sysname"] == "Windows"
}
