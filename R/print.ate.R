### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: okt 29 2019 (13:36) 
##           By: Brice Ozenne
##     Update #: 346
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.ate (documentation)
#' @title Print Average Treatment Effects
#' @description Print average treatment effects.
#' @name print.ate
#' 
#' @param x object obtained with function \code{ate}
#' @param ... passed to summary
#'
#' @seealso
#' \code{\link{summary.ate}} to obtained a more detailed output
#' \code{\link{confint.ate}} to compute confidence intervals/bands.
#' \code{\link{ate}} to compute the average treatment effects.

## * print.ate (code)
#' @rdname print.ate
#' @method print ate
#' @export
print.ate <- function(x, ...){
    summary(x, short = TRUE, ...)
}


#----------------------------------------------------------------------
### print.ate.R ends here
