# .onLoad()/.onUnload() is used when there is a NAMESPACE file
.onLoad <- function(libname, pkgname) {
  library.dynam(pkgname, package=pkgname, lib.loc=.libPaths())
}

.onUnload <- function(libpath) {
  pkgname <- "HaarSeg";
  library.dynam.unload(chname=pkgname, libpath=libpath)
}





