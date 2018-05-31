### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: maj 31 2018 (11:53) 
##           By: Brice Ozenne
##     Update #: 47
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
#' @title Turn ate object into a \code{data.table}
#' @description Turn ate object into a \code{data.table}
#' @param x object obtained with function \code{ate}
#' @param keep.rownames not used
#' @param se should standard errors/quantile for confidence bands be displayed?
#' @param ... not used
#' @export
as.data.table.ate <- function(x, keep.rownames = FALSE, se = TRUE, ...){

    ## ** ate
    out1 <- data.table(type = "ate",
                       level = x$meanRisk[[1]],
                       time = x$meanRisk[["time"]],
                       value = x$meanRisk[["meanRisk"]])
    if(x$se && !is.null(x$conf.level)){
        out1$se <- x$meanRisk[["se"]]
        out1$lower <- x$meanRisk[["lower"]]
        out1$upper <- x$meanRisk[["upper"]]
    }
    if(x$band && !is.null(x$conf.level)){
        out1$lowerBand <- x$meanRisk[["lowerBand"]]
        out1$upperBand <- x$meanRisk[["upperBand"]]
    }
    
    ## ** diff ate
    out2 <- data.table(type = "diffAte",
                       level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                       time = x$riskComparison[["time"]],
                       value = x$riskComparison[["diff"]])
    if(x$se && !is.null(x$conf.level)){
        if(se){out2$se <- x$riskComparison[["diff.se"]]}
        out2$lower <- x$riskComparison[["diff.lower"]]
        out2$upper <- x$riskComparison[["diff.upper"]]
    }
    if(x$band && !is.null(x$diffAte.transform)){
        if(se){out2$lowerBand <- x$riskComparison[["diff.lowerBand"]]}
        out2$upperBand <- x$riskComparison[["diff.upperBand"]]
    }

    ## ** ratio ate
    out3 <- data.table(type = "ratioAte",
                       level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                       time = x$riskComparison[["time"]],
                       value = x$riskComparison[["ratio"]])
    
    if(x$se && !is.null(x$ratioAte.transform)){
        if(se){out3$se <- x$riskComparison[["ratio.se"]]}
        out3$lower <- x$riskComparison[["ratio.lower"]]
        out3$upper <- x$riskComparison[["ratio.upper"]]
    }
    if(x$band && !is.null(x$ratioAte.transform)){
        if(se){out3$lowerBand <- x$riskComparison[["ratio.lowerBand"]]}
        out3$upperBand <- x$riskComparison[["ratio.upperBand"]]
    }
    
    ## export
    return(rbind(out1,out2,out3))
  
}



######################################################################
### as.data.table.ate.R ends here
