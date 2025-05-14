### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 21 2020 (11:10) 
## Version: 
## Last-Updated: sep 11 2024 (18:31) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @export
nobs.CauseSpecificCox <- function(object,...){
    return(NROW(object$response))
}

##----------------------------------------------------------------------
### nobs.R ends here
