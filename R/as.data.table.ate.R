### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: okt 11 2019 (13:45) 
##           By: Brice Ozenne
##     Update #: 105
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
#' @param estimator [character] The type of estimator relative to which the estimates should be output. 
#' @param ... Not used.
#'

## * as.data.table.ate (code)
#' @rdname as.data.table.ate
#' @export
as.data.table.ate <- function(x, keep.rownames = FALSE, se = TRUE, estimator = x$estimator, ...){

    estimator <- match.arg(estimator, choices =  x$estimator, several.ok = TRUE)
    landmark <- all(attr(x$estimator,"full")=="GformulaTD")
    col.time <- switch(as.character(landmark),
                       "TRUE" = c("time","landmark"),
                       "FALSE" = "time")
    out <- NULL
    
    for(iE in 1:length(estimator)){ ## iE <- 1
        iEstimator <- estimator[iE]

        ## ** which columns to keep
        keep.col <- NULL
        if(x$se[[1]] && !is.null(x$conf.level)){
            if(se){
                keep.col <- c(keep.col, "se")
            }
            keep.col <- c(keep.col, "lower", "upper", "p.value")
            x$meanRisk[, paste0("meanRisk.",iEstimator,".p.value") := as.numeric(NA)]
        }
        if(x$band[[1]] && !is.null(x$conf.level)){
            if(se){
                keep.col <- c(keep.col, "quantileBand")
            }
            keep.col <- c(keep.col, "lowerBand", "upperBand")
        }
        
    
        ## ** ate
        out1 <- cbind(type = "ate",
                      level = x$meanRisk[[1]],
                      x$meanRisk[,.SD,.SDcols = col.time],
                      value = x$meanRisk[[paste0("meanRisk.",iEstimator)]])
        if(!is.null(x$boot)){
            out1[,c("value.boot") := x$meanRisk[[paste0("meanRisk.",iEstimator,".boot")]]]
        }
        if(!is.null(keep.col)){
            out1[, c(keep.col) := x$meanRisk[, .SD, .SDcols = paste0("meanRisk.",iEstimator,".",keep.col)]]
        }
        out1[, c("estimator") := iEstimator]
        
        ## ** diff ate
        out2 <- cbind(type = "diffAte",
                      level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                      x$riskComparison[,.SD,.SDcols = col.time],
                      value = x$riskComparison[[paste0("diff.",iEstimator)]])
        if(!is.null(x$boot)){
            out2[,c("value.boot") := x$riskComparison[[paste0("diff.",iEstimator,".boot")]]]
        }
        if(!is.null(keep.col)){
            out2[, c(keep.col) := x$riskComparison[, .SD, .SDcols = paste0("diff.",iEstimator,".",keep.col)]]
        }
        out2[, c("estimator") := iEstimator]
    
        ## ** ratio ate
        out3 <- cbind(type = "ratioAte",
                      level = paste0(x$riskComparison[[1]],".",x$riskComparison[[2]]),
                      x$riskComparison[,.SD,.SDcols = col.time],
                      value = x$riskComparison[[paste0("ratio.",iEstimator)]])
        if(!is.null(x$boot)){
            out3[,c("value.boot") := x$riskComparison[[paste0("ratio.",iEstimator,".boot")]]]
        }
        if(!is.null(keep.col)){
            out3[, c(keep.col) := x$riskComparison[, .SD, .SDcols = paste0("ratio.",iEstimator,".",keep.col)]]
        }
        out3[, c("estimator") := iEstimator]

        ## ** assemble
        out <- rbind(out,rbind(out1,out2,out3))
    }
    
    ## export
    return(out)
  
}



######################################################################
### as.data.table.ate.R ends here
