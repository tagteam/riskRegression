### print.ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (15:28) 
## Version: 
## Last-Updated: jul 18 2018 (13:50) 
##           By: Brice Ozenne
##     Update #: 15
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
#' @param efficient [logical] should the standard errors be corrected to improve efficiency
#' in presence of censoring.
#' @param nuisance.iid [logical] should the standard errors
#' accounting for the uncertainty in the outcome and propensity score models be output?
#' @param ... Passed to print.
#' 


## * print.ateRobust (code)
#' @rdname print.ateRobust
#' @method print ateRobust
#' @export
print.ateRobust <- function(x, digits = 3, efficient = x$efficient, nuisance.iid = x$nuisance.iid,...){

    if(nuisance.iid && (x$nuisance.iid == FALSE)){
        stop("Argument \'nuisance.iid\' must be FALSE when ateRobust has been called using \'iid\'=FALSE \n")
    }

    ## columns to keep
    if(efficient && nuisance.iid){
        keep.cols <- c("Gformula2","IPWefficient2","AIPWefficient2")
    }else if(efficient){
        keep.cols <- c("Gformula","IPWefficient","AIPWefficient")
    }else if(nuisance.iid){
        keep.cols <- c("Gformula2","IPWnaive2","AIPWnaive2")
    }else{
        keep.cols <- c("Gformula","IPWnaive","AIPWnaive")
    }

    ## extract
    x <- list(ate.value = x$ate.value[,keep.cols],
              ate.se = x$ate.se[,keep.cols])

    ## output
    print(x, digits=digits, ...)
    invisible(x)
}



######################################################################
### print.ateRobust.R ends here
