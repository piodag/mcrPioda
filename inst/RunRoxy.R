##############################################################################
##
## RunRoxy.R
##
## Generate Roxygen documentation from source code.
## Note that Roxygen has problems with S4 classes and therefore all
## class definition files (MCResult*.R) have to be removed from R directory.
##
###############################################################################

library(roxygen2)
roxygenize(package.dir="mcrPioda")
