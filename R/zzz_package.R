# This is going to run when the user loads the package
#
.onAttach <- function(libname, pkgname){

  package_ver <- utils::packageDescription("chouca",
                                           fields = "Version")

  devel_message <- ""
  if ( grepl("99$", package_ver, perl = TRUE) ) {
    devel_message <- " (development version, use at your own risk!)"
  }

  packageStartupMessage({
    paste0("This is chouca ", package_ver, devel_message, "\n",
           "See vignette(\"chouca-package\") for a user guide or ?chouca for an overview")
    }, appendLF = TRUE)

}
