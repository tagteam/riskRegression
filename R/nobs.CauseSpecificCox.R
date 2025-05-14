### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 21 2020 (11:10) 
## Version: 
## Last-Updated: May 14 2025 (08:52) 
##           By: Thomas Alexander Gerds
##     Update #: 16
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @method nobs CauseSpecificCox 
##' @export
nobs.CauseSpecificCox <- function(object,...){
    return(NROW(object$response))
}
#' @export
nobs.coxph <- function(object,...){
    return(object$n)
}
#' @export
nobs.phreg <- function(object,...){
    return(NROW(object$time))
}


##----------------------------------------------------------------------
### nobs.R ends here
