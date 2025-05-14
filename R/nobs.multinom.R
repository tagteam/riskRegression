### nobs.multinom.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 14 2025 (08:52) 
## Version: 
## Last-Updated: May 14 2025 (08:52) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @method nobs multinom 
##' @export
nobs.multinom <- function(object,...){
    NROW(object$residuals)
}
######################################################################
### nobs.multinom.R ends here
