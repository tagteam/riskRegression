### terms.phreg.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:40) 
## Version: 
## Last-Updated: Apr 27 2025 (07:40) 
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
## * terms
## ** terms.phreg
#' @title Extract terms for phreg objects
#' @description Extract terms for phreg objects
#' @param x a phreg object.
#' @param ... not used.
#' 
#' @method terms phreg
terms.phreg <- function(x, ...){
    stats::terms(x$formula)
}
######################################################################
### terms.phreg.R ends here
