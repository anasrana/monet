.onLoad <- function(libname, pkgname) {
  if (length(stanmodels) != 0) {
    modules <- paste0("stan_fit4", names(stanmodels), "_mod")
    for (m in modules) loadModule(m, what = TRUE)
  } else {
    message("No stan programs to compile were found.")
  }
}

.onAttach <- function(libname, pkgname) {
  if (!interactive()) {
    return()
  }

  v = packageVersion("monet")
  d = read.dcf(system.file("DESCRIPTION", package="monet"),
               fields = c("Packaged", "Built", "Revision"))
  startupMessage = ""
    packageStartupMessage(
      "Package monet loaded"
    )
}
