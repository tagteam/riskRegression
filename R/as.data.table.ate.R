### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: jun 13 2018 (14:57) 
##           By: Brice Ozenne
##     Update #: 73
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.table.ate (documentation)
#' @title Turn ate Object Into a \code{data.table}
#' @description Turn ate object into a \code{data.table}.
#' @name as.data.table.ate
#' 
#' @param x object obtained with function \code{ate}
#' @param keep.rownames Not used.
#' @param se [logical] Should standard errors/quantile for confidence bands be displayed?
#' @param ... Not used.
#'

## * as.data.table.ate (code)
#' @rdname as.data.table.ate
#' @export
as.data.table.ate <- function(x, keep.rownames = FALSE, se = TRUE, ...){

    ## ** which columns to keep
    keep.col <- NULL
    if(x$se && !is.null(x$conf.level)){
        if(se){
            keep.col <- c(keep.col, "se")
        }
        keep.col <- c(keep.col, "lower", "upper")
    }
    if(x$band && !is.null(x$conf.level)){
        if(se){
            keep.col <- c(keep.col, "quantileBand")
        }
        keep.col <- c(keep.col, "lowerBand", "upperBand")
    }
    
    keep.cols <- c("lower","upper")
    
    ## ** ate
    out1 <- data.table(type = "ate",
                       level = x$meanRisk[[1]],
                       time = x$meanRisk[["time"]],
                       value = x$meanRisk[["meanRisk"]])
    if(!is.null(x$boot)){
        out1[,c("value.boot") := x$meanRisk$meanRisk.boot]
    }
    out1[, c(keep.col) := x$meanRisk[, .SD, .SDcols = paste0("meanRisk.",keep.col)]]
    
    ## ** diff ate
    out2 <- data.table(type = "diffAte",
                       level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                       time = x$riskComparison[["time"]],
                       value = x$riskComparison[["diff"]])
    if(!is.null(x$boot)){
        out2[,c("value.boot") := x$riskComparison$diff.boot]
    }

    out2[, c(keep.col) := x$riskComparison[, .SD, .SDcols = paste0("diff.",keep.col)]]

    ## ** ratio ate
    out3 <- data.table(type = "ratioAte",
                       level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                       time = x$riskComparison[["time"]],
                       value = x$riskComparison[["ratio"]])
    
    out3[, c(keep.col) := x$riskComparison[, .SD, .SDcols = paste0("ratio.",keep.col)]]
    if(!is.null(x$boot)){
        out3[,c("value.boot") := x$riskComparison$ratio.boot]
    }

    ## export
    return(rbind(out1,out2,out3))
  
}



######################################################################
### as.data.table.ate.R ends here
