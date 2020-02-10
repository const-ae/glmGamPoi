

is_macos <- function(){
  Sys.info()["sysname"] == "Darwin"
}
