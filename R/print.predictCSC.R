### print.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 11 2017 (10:01) 
## Version: 
## last-updated: Sep  4 2017 (17:31) 
##           By: Thomas Alexander Gerds
##     Update #: 61
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Print predictions from a Cause-specific Cox proportional hazard regression
#' @description Print predictions from a Cause-specific Cox proportional hazard regression
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @param digits integer indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @examples
#' ## no strata
#' d <- sampleData(1e2, outcome = "competing.risks")
#' m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)
#' pred.CSC <- predict(m.CSC, time = 1:5, cause = 1,
#'                       se = TRUE, keep.newdata = TRUE)
#'
#' pred.CSC
#' print(pred.CSC, ci = TRUE)
#'
#' ## strata
#' library(survival)
#' m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6,
#'               data = d)
#' pred.SCSC <- predict(m.SCSC, time = 1:5, cause = 1,
#'                se = TRUE, keep.newdata = TRUE, keep.strata = TRUE)
#' pred.SCSC
#' print(pred.SCSC, ci = TRUE)
#' 
#' @method print predictCSC
#' @export
print.predictCSC <- function(x,
                             digits = 3, ...){
        out <- as.data.table(x)
        print(out,digits=digits,...)
        invisible(out)
}


## `[.predictCSC` <- function(x, i, j, drop = FALSE){

##     if(missing(i)){
##         i <- 1:NROW(x$absRisk)
##     }
##     if(missing(j)){
##         j <- 1:NCOL(x$absRisk)
##     }

##     return(x$absRisk[i,j,drop = drop])
## }

#----------------------------------------------------------------------
### print.predictCSC.R ends here
