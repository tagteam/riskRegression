### print.ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (15:28) 
## Version: 
## Last-Updated: Sep  7 2018 (11:12) 
##           By: Thomas Alexander Gerds
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.ateRobust (documentation)
#' @title Print Average Treatment Effect
#' @description Print average treatment effect.
#' @name print.ateRobust
#' 
#' @param x object obtained with the function \code{ateRobust}.
#' @param digits [integer, >0] indicating the number of decimal places.
#' @param augment.cens [logical] should the standard errors account for the augmentation term
#' in direction of the model for the censoring mechanism.
#' accounting for the uncertainty in the outcome and propensity score models be output
#' @param ... Passed to print.
#' 


## * print.ateRobust (code)
#' @rdname print.ateRobust
#' @method print ateRobust
#' @export
print.ateRobust <- function(x, digits = 3, augment.cens = x$augment.cens,...){


    ## columns to keep
    ## keep.cols <- c("Gformula","IPTW.AIPCW","AIPTW.AIPCW")

    ## extract
    x <- list(ate.value = x$ate.value,
              ate.se = x$ate.se)

    ## output
    print(x, digits=digits, ...)
    invisible(x)
}



######################################################################
### print.ateRobust.R ends here
