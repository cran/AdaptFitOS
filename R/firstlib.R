.First.lib <- function(lib,pkg)
{
   library.dynam("AdaptFitOS",pkg,lib)
   ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
   ver <- as.character(ver)
#   cat("AdaptFitOS",ver,"loaded\n")
#   cat("Note: This package replaces functions of package ConfBands which will lead to conflicts when both packages are used at the same time.\n")
  hello <- paste("This is AdaptFitOS ",ver,". Type 'help(\"AdaptFitOS-package\")' for an overview.",sep="")
  packageStartupMessage(hello)
}

