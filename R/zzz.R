# global reference to scipy (will be initialized in .onLoad)
quaternion <- NULL
numpy <- NULL
math <- NULL



#' @importFrom reticulate import
.onLoad <- function(libname, pkgname) {

  # user_permission <- utils::askYesNo("Install miniconda? downloads 50MB and takes time")
  #
  # if (isTRUE(user_permission)) {
  #   reticulate::install_miniconda()
  # } else {
  #   message("You should run `reticulate::install_miniconda()` before using this package")
  # }

  # use superassignment to update global reference to scipy
  numpy <<- reticulate::import("numpy", delay_load = TRUE)
  math <<- reticulate::import("math", delay_load = TRUE)
  quaternion <<- reticulate::import("quaternion", delay_load = TRUE)

  # Global
  #deg2rad <- pi / 180
  #reticulate::source_python("R/quaternions.py", convert = FALSE)

}
