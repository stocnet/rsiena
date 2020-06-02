#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: zzz.r
# *
# * Description: This module contains the code for package attachment and
# * detachment.
# *****************************************************************************/

#dllpath <- ''
##@imagepath Objects/Path Path for image of Siena
imagepath <- ""

##csvpath <- ''
## or .onAttach?
##@pkgpath Objects/Path Path for installation of siena.exe
pkgpath <- ""

##@pkgname Objects/Package name of package (RSiena or RSienaTest)
pkgname <- ""

##@pkgversion Objects/Package version of package
pkgvers <- ""

##@.onLoad Miscellaneous Start-up processing
.onLoad <- function(libname, pkgname) {
#   dllpath <<- if (nzchar(.Platform$r_arch))
##       file.path(libname, pkgname, "libs", .Platform$r_arch,
#                 paste('RSiena', .Platform$dynlib.ext, sep=''))
#   else
#       file.path(libname,pkgname, "libs",
 #                paste('RSiena', .Platform$dynlib.ext, sep=''))
#  #  data('sysdata',package='RSiena')
  imagepath <<- file.path(libname, pkgname, paste("ilcampo.gif"))
  pkgpath <<- file.path(libname, pkgname)
  pkgname <<- pkgname
  pkgvers <<- utils::packageDescription(pkgname, fields = c("Version", "Date"))
#cat(pkgname,pkgpath,'\n')
 ##  csvpath<<- file.path(libname,pkgname)
 #   library.dynam("RSiena",package=pkgname)
 #  cat (libname,pkgname,'\n')
}

##@.onUnload Miscellaneous Unload processing
.onUnload <- function(libpath) {
  library.dynam.unload(pkgname, libpath)
}

#.Last.lib <- function(libpath)
#{
#    cat(libpath,'\n')
#}


