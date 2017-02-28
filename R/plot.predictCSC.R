### plot.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 27 2017 (10:47) 
## Version: 
## last-updated: Feb 28 2017 (18:47) 
##           By: Thomas Alexander Gerds
##     Update #: 34
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Plot predictions from a Cause-specific Cox proportional hazard regression
#' @description Plot predictions from a Cause-specific Cox proportional hazard regression
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param groupBy The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot Logical. Should the graphic be plotted.
#' @param conf.level confidence level of the interval.
#' @param digit integer indicating the number of decimal places
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @examples
#' ## no strata
#' d <- sampleData(1e2, outcome = "competing.risks")
#' m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)
#' 
#' pred.CSC <- predict(m.CSC, time = 1:5, cause = 1)
#' plot(pred.CSC)
#' 
#' pred.CSC <- predict(m.CSC, newdata = d[1:3,],
#'                     time = 1:5, cause = 1, se = TRUE, keep.newdata = TRUE)
#' plot(pred.CSC, groupBy = "covariates")
#'
#' ## strata
#' m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6,
#' data = d)
#' pred.SCSC <- predict(m.SCSC, time = 1:3,
#' cause = 1, se = TRUE, keep.newdata = TRUE, keep.strata = TRUE)
#' plot(pred.SCSC)
#' plot(pred.SCSC, groupBy = "strata")
#'
#' @method plot predictCSC
#' 
#' @export
plot.predictCSC <- function(x,
                            ci = FALSE,
                            groupBy = "row",
                            reduce.data = FALSE,
                            plot = TRUE,
                            conf.level = 0.95,
                            digit = 2, ...){
    
    ## initialize and check        
    possibleGroupBy <- c("row","covariates","strata")
    if(groupBy %in% possibleGroupBy == FALSE){
        stop("argument \"groupBy\" must be in \"",paste(possibleGroupBy, collapse = "\" \""),"\"\n")
    }
    
    if(groupBy == "covariates" && ("newdata" %in% names(x) == FALSE)){
        stop("argument \'groupBy\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling predictCox \n")
    }
    if(groupBy == "strata" && ("strata" %in% names(x) == FALSE)){
        stop("argument \'groupBy\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling predictCox \n")
    }

    if(ci && "absRisk.se" %in% names(x) == FALSE){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }

    ## display
    newdata <- copy(x$newdata)
    if(!is.null(newdata) && reduce.data){
        test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
        if(any(test)){
            newdata[, (names(test)[test]):=NULL]
        }        
    }

    gg.res <- predict2plot(outcome = x$absRisk,
                           outcome.se = if(ci){x$absRisk.se}else{NULL},
                           newdata = newdata,
                           strata = x$strata,
                           times = x$times,
                           digit = digit,
                           name.outcome = "absoluteRisk", # must not contain space to avoid error in ggplot2
                           conf.level = conf.level,
                           groupBy = groupBy)

    gg.res <- gg.res + xlab("absolute risk")
    
    if(plot){
        print(gg.res)
    }
    
    return(invisible(list(plot = gg.res,
                          data = newdata)))
}

#----------------------------------------------------------------------
### plot.predictCSC.R ends here
