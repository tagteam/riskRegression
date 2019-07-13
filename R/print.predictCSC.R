### print.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 11 2017 (10:01) 
## Version: 
## last-updated: Jul 13 2019 (10:54) 
##           By: Thomas Alexander Gerds
##     Update #: 77
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.predictCSC (documentation)
#' @title Print Predictions From a Cause-specific Cox Proportional Hazard Regression
#' @description Print predictions from a Cause-specific Cox proportional hazard regression.
#' @name print.predictCSC
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @param digits [integer, >0] indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @details to display confidence intervals/bands,
#' the \code{confint} method needs to be applied on the object.
#'
#' @seealso
#' \code{\link{confint.predictCSC}} to compute confidence intervals/bands.
#' \code{\link{predict.CauseSpecificCox}} to compute the predicted risks.

## * print.predictCSC (code)
#' @rdname print.predictCSC
#' @method print predictCSC
#' @export
print.predictCSC <- function(x, digits = 3, ...){
        x <- as.data.table(x)
        print(x,digits=digits,...)
        invisible(x)
}

#----------------------------------------------------------------------
### print.predictCSC.R ends here
