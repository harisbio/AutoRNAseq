#' @importFrom shiny addResourcePath

.onLoad <- function(libname, pkgname) {
  # Create link to logo
  options(shiny.maxRequestSize = 8192 * 1024^2)
  shiny::addResourcePath("AutoRNAseq", system.file("www", package = "AutoRNAseq"))

}
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription("AutoRNAseq", fields = "Version")
  msg <- paste0("Welcome to AutoRNAseq v", pkgVersion, "\n")
  packageStartupMessage(msg)

}

