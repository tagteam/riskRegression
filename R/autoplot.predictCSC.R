### autoplot.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 27 2017 (10:47) 
## Version: 
## last-updated: okt  3 2017 (17:19) 
##           By: Brice Ozenne
##     Update #: 65
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
#' @param object object obtained with the function \code{predictCox}.
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param band Logical. If \code{TRUE} display the confidence bands for the predictions.
#' @param group.by The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}. 
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot Logical. Should the graphic be plotted.
#' @param digits integer indicating the number of decimal places
#' @param alpha transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @examples
#' ## no strata
#' d <- sampleData(1e2, outcome = "competing.risks")
#' m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)
#' 
#' pred.CSC <- predict(m.CSC, newdata = d[1:2,], time = 1:5, cause = 1)
#' autoplot(pred.CSC)
#' 
#' pred.CSC <- predict(m.CSC, newdata = d[1:3,],
#'                     time = 1:5, cause = 1, se = TRUE, keep.newdata = TRUE)
#' 
#'
#' ## strata
#' library(survival)
#' m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6,
#' data = d)
#' pred.SCSC <- predict(m.SCSC, time = 1:3, newdata = d[1:4,],
#' cause = 1, se = TRUE, keep.newdata = TRUE, keep.strata = TRUE)
#' autoplot(pred.SCSC, group.by = "strata")
#'
#' @method autoplot predictCSC
#' 
#' @export
autoplot.predictCSC <- function(object,
                            ci = FALSE,
                            band = FALSE,
                            group.by = "row",
                            reduce.data = FALSE,
                            plot = TRUE,
                            digits = 2, alpha = NA, ...){
  
    ## initialize and check        
	group.by <- match.arg(group.by, c("row","covariates","strata"))
  
    if(group.by == "covariates" && ("newdata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling predictCox \n")
    }
    if(group.by == "strata" && ("strata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling predictCox \n")
    }
  
    if(ci && object$se == FALSE){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }
    if(band && object$band == FALSE){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set argment \'nsim.band\' to a positive integer when calling predict.CauseSpecificCox \n")
    }
    
  ## display
  newdata <- copy(object$newdata)
  if(!is.null(newdata) && reduce.data){
    test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
    if(any(test)){
      newdata[, (names(test)[test]):=NULL]
    }        
  }

    dataL <- predict2melt(outcome = object$absRisk, ci = ci, band = band, 
                          outcome.lower = if(ci){object$absRisk.lower}else{NULL},
                          outcome.upper = if(ci){object$absRisk.upper}else{NULL},
                          outcome.lowerBand = if(band){object$absRisk.lowerBand}else{NULL},
                          outcome.upperBand = if(band){object$absRisk.upperBand}else{NULL},
                          newdata = newdata,
                          strata = object$strata,
                          times = object$times,
                          name.outcome = "absoluteRisk",
                          group.by = group.by,
                          digits = digits
                          )

    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = "absoluteRisk", # must not contain space to avoid error in ggplot2
                           ci = ci,
                           band = band,
                           group.by = group.by,
                           conf.level = object$conf.level,
                           alpha = alpha,
                           origin = min(object$time)
                           )
    
  gg.res$plot <- gg.res$plot + ggplot2::ylab("absolute risk")
  
  if(plot){
    print(gg.res$plot)
  }
  
  return(invisible(gg.res))
}

#----------------------------------------------------------------------
### autoplot.predictCSC.R ends here
