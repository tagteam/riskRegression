### nobs.R ---
## ----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 21 2020 (11:10)
## Version:
## Last-Updated: apr 21 2020 (11:16)
##           By: Brice Ozenne
##     Update #: 11
## ----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
## ----------------------------------------------------------------------
##
### Code:

nobs.CauseSpecificCox <- function(object, ...) {
  return(NROW(object$response))
}
nobs.coxph <- function(object, ...) {
  return(object$n)
}
nobs.phreg <- function(object, ...) {
  return(NROW(object$time))
}

## ----------------------------------------------------------------------
### nobs.R ends here
