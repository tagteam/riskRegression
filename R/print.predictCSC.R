### print.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 11 2017 (10:01) 
## Version: 
## last-updated: feb 27 2017 (11:47) 
##           By: Brice Ozenne
##     Update #: 57
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
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param conf.level confidence level of the interval.
#' @param digit integer indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @examples
#' ## no strata
#' d <- sampleData(1e2, outcome = "competing.risks")
#' m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)
#' pred.CSC <- predict(m.CSC, time = 1:5, cause = 1, se = TRUE, keep.newdata = TRUE)
#'
#' pred.CSC
#' print(pred.CSC, ci = TRUE)
#'
#' ## strata
#' m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6, data = d)
#' pred.SCSC <- predict(m.SCSC, time = 1:5, cause = 1, se = TRUE, keep.newdata = TRUE, keep.strata = TRUE)
#' pred.SCSC
#' print(pred.SCSC, ci = TRUE)
#' 
#' @method print predictCSC
#' @export
print.predictCSC <- function(x,
                             ci = FALSE,
                             reduce.data = FALSE,
                             conf.level = 0.95,
                             digit = 3, ...){

    ## check the presence of the standard errors
    if(ci && "absRisk.se" %in% names(x)==FALSE){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }

    ## remove covariates that have the same value for all observations
    newdata <- copy(x$newdata)
    if(!is.null(newdata) && reduce.data){
        test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
        if(any(test)){
            newdata[, (names(test)[test]):=NULL]
        }        
    }

    ##
    ls.print <- predict2print(x$absRisk,
                              outcome.se = if(ci){x$absRisk.se}else{NULL},
                              newdata = x$newdata,
                              strata = x$strata,
                              times = x$times,
                              digit = digit,
                              name.outcome = "absolute risk",
                              conf.level = conf.level,
                              lower = 0, upper = Inf)

    resDisplay <- ls.print$resDisplay
    res <- ls.print$res

    ## display
    cat("absolute risk (columns: ",
        if(!is.null(x$strata)){"strata|"},
        if(!is.null(newdata)){"covariate|"},
        "time, rows: individuals",
        if(ci){paste0(", [.;.] ",100*conf.level,"% confidence interval) \n")}else{") \n"},
        sep = "")
    print(noquote(resDisplay), ...)

    ## export
    return(invisible(res))

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
