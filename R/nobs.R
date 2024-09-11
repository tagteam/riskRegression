### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 21 2020 (11:10) 
## Version: 
## Last-Updated: Jun 13 2024 (19:54) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @export
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
