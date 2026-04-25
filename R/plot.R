### plot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr 23 2026 (19:24) 
## Version: 
## Last-Updated: Apr 25 2026 (12:46) 
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

## * plot.predictCox (code)
##' @describeIn autoplot.predictCox Graphical Display of Predictions From a Cox Model
##' @export
plot.predictCox <- function(x, ...){
    out <- autoplot.predictCox(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object (e.g. when using qqtest)
        print(out$plot)
    }
    return(invisible(out))

}

## * plot.predictCSC (code)
##' @describeIn autoplot.predictCSC Graphical Display of Predictions From a Cause-Specific Cox Model
##' @export
plot.predictCSC <- function(x, ...){
    out <- autoplot.predictCSC(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object (e.g. when using qqtest)
        print(out$plot)
    }
    return(invisible(out))

}

## * plot.ate (code)
##' @describeIn autoplot.ate Graphical Display of Average Risks
##' @export
plot.ate <- function(x, ...){
    out <- autoplot.ate(x, ...)
    if(!is.null(out$plot)){
        ## if ggplot then the graph is stored in out and should be displayed here
        ## otherwise the graph has already been displayed but could not be stored in the object (e.g. when using qqtest)
        print(out$plot)
    }
    return(invisible(out))
}

##----------------------------------------------------------------------
### plot.R ends here
