
.onLoad <- function(libname, pkgname) {
ver<- utils::packageDescription(pkgname, libname, fields=c("Version", "Date"))
  hello <- paste("This is AdaptFitOS ",ver[1],". Type 'help(\"AdaptFitOS-package\")' for an overview.",sep="")
  packageStartupMessage(hello)
}

