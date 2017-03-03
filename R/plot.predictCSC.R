### plot.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 27 2017 (10:47) 
## Version: 
## last-updated: Mar  3 2017 (11:41) 
##           By: Thomas Alexander Gerds
##     Update #: 37
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
#' @param digits integer indicating the number of decimal places
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
                            digits = 2, ...){
    plotframe <- as.data.table(x)
}

#----------------------------------------------------------------------
### plot.predictCSC.R ends here
