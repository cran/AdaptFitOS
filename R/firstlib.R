.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("AdaptFitOS ",ver," loaded.\n Type 'help(\"AdaptFitOS-package\")' for an overview.\n", domain = NULL,  appendLF = TRUE)
}