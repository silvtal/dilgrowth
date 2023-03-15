.onUnload <- function (libpath) {
  ## Whenever you use C++ code in your package, you need to clean up after
  ## yourself when your package is unloaded
  library.dynam.unload("dilgrowth", libpath)
}
