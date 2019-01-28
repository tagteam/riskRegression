### print.influenceTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  1 2018 (13:35) 
## Version: 
## Last-Updated: Jan 28 2019 (15:51) 
##           By: Thomas Alexander Gerds
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.influenceTest (documentation)
#' @title Output of the DIfference Between Two Estimates
#' @description Output of the difference between two estimates.
#' @rdname print.influenceTest
#' 
#' @param x object obtained with the function \code{influenceTest}.
#' @param digits [integer, >0] indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @details to display confidence intervals/bands,
#' the \code{confint} method needs to be applied on the object.
#'
#' @seealso
#' \code{\link{confint.influenceTest}} to compute confidence intervals/bands.
#' \code{\link{influenceTest}} to perform the comparison.

## * print.predictCSC (code)
#' @rdname print.influenceTest
#' @method print influenceTest
#' @export
print.influenceTest <- function(x, digits = 3, ...){

    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table
    
    cat("        Comparison of two estimates of the ",x$type," at time ",paste0(round(x$time, digits), collapse = " "),"\n\n", sep = "")

    cat(x$name[1],":\n",sep="")
    print(x$call[[1]])
    cat(x$name[2],":\n",sep="")
    print(x$call[[2]])

    dt.tempo <- as.data.table(x)
    order.col <- setdiff(names(dt.tempo),c("lower","upper","p.value","quantileBand","lowerBand","upperBand"))

    ## round and merge column containing CI and CB
    test.numeric <- unlist(lapply(dt.tempo, is.numeric))
    numeric.col <- setdiff(names(dt.tempo)[test.numeric],"p.value")
    ## dt.tempo[, c(numeric.col) := round(.SD, digits = digits) , .SDcols = numeric.col]

    ## x$se
    if(!is.null(x$conf.level)){
        dt.tempo[, c("conf.interval") := paste0("[",lower," ; ",upper,"]")]
        dt.tempo[,c("lower","upper") := NULL]
        order.col <- c(order.col,"conf.interval","p.value")
    }
    if(x$band && !is.null(x$conf.level)){
        dt.tempo[, c("conf.band") := paste0("[",lowerBand," ; ",upperBand,"]")]
        dt.tempo[,c("lowerBand","upperBand") := NULL]
        order.col <- c(order.col,"quantileBand","conf.band")
    }
    ## round p.value
    if(!is.null(x$conf.level)){
        dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
    }

    ## print
    data.table::setcolorder(dt.tempo, neworder = order.col)
    print(dt.tempo,digits=digits,...)

    ## export
    return(invisible(x))
}


######################################################################
### print.influenceTest.R ends here
