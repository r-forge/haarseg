# .onLoad()/.onUnload() is used when there is a NAMESPACE file
.onLoad <- function(libname, pkgname) {
  library.dynam(pkgname, package=pkgname)
}

.onUnload <- function(libpath) {
  library.dynam.unload(libpath=libpath)
}

# .onLoad()/.onUnload() is used when there is *no* NAMESPACE file
# .First.lib <- function(libname, pkgname) {
#   library.dynam(pkgname, package=pkgname)
# }
# 
# .Last.lib <- function(libpath) {
#   library.dynam.unload(libpath=libpath)
# }





