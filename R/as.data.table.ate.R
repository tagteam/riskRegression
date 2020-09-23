### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: sep 23 2020 (14:00) 
##           By: Brice Ozenne
##     Update #: 157
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
#' @param se [logical] Should the standard errors be output?
#' @param ci [logical] Should the confidence intervals be output?
#' @param band [logical] Should the confidence bands be output?
#' @param quantileBand [logical] Should the quantiles for confidence bands be output?
#' @param p.value [logical] Should the p-values/adjusted p-values be output?
#' @param estimator [character] The type of estimator relative to which the estimates should be output. 
#' @param ... Not used.
#'

## * as.data.table.ate (code)
#' @rdname as.data.table.ate
#' @export
as.data.table.ate <- function(x, estimator = x$estimator,
                              se = x$se, ci = x$ci, band = x$band, quantileBand = x$band, p.value = x$p.value,
                              keep.rownames = FALSE, ...){

    estimator <- match.arg(estimator, choices =  x$estimator, several.ok = TRUE)

    landmark <- all(attr(x$estimator,"full")=="GFORMULATD")
    col.time <- switch(as.character(landmark),
                       "TRUE" = c("time","landmark"),
                       "FALSE" = "time")
    out <- NULL
    allContrasts <- x$allContrasts
    contrasts <- attr(allContrasts,"contrasts")
    
    meanRisk <- data.table::copy(x$meanRisk)
    diffRisk <- data.table::copy(x$diffRisk)
    ratioRisk <- data.table::copy(x$ratioRisk)

    ## ** meanRisk
    out1 <- cbind(type = "meanRisk",
                  estimator = x$meanRisk$estimator,
                  time = x$meanRisk$time,
                  level = x$meanRisk$treatment,
                  x$meanRisk[,.SD,.SDcols = setdiff(names(x$meanRisk),c("estimator","time","treatment"))])
    if(p.value){
        out1$p.value <- as.numeric(NA)
    }


    ## ** diffRisk
    out2 <- cbind(type = "diffRisk",
                  estimator = x$diffRisk$estimator,
                  time = x$diffRisk$time,
                  level = paste0(x$diffRisk$A,".",x$diffRisk$B),
                  x$diffRisk[,.SD,.SDcols = setdiff(names(x$diffRisk),c("estimator","time","A","B","estimate.A","estimate.B"))])

    ## ** ratioRisk
    out3 <- cbind(type = "ratioRisk",
                  estimator = x$ratioRisk$estimator,
                  time = x$ratioRisk$time,
                  level = paste0(x$ratioRisk$A,".",x$ratioRisk$B),
                  x$ratioRisk[,.SD,.SDcols = setdiff(names(x$ratioRisk),c("estimator","time","A","B","estimate.A","estimate.B"))])
        
    
    ## export
    return(rbind(out1,out2,out3))
  
}



######################################################################
### as.data.table.ate.R ends here
