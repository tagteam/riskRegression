### print.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 15 2017 (17:36) 
## Version: 
## last-updated: Mar  3 2017 (09:55) 
##           By: Brice Ozenne
##     Update #: 143
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Print predictions from a Cox model
#' @description Print predictions from a Cox model
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @inheritParams predictCox
#' @param digits integer indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @examples
#' library(survival)
#' library(rms)
#' 
#' set.seed(10)
#' d <- sampleData(1e2, outcome = "survival")
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                data = d, x = TRUE, y = TRUE)
#' predictCox(m.cox)
#'
#' pred <- predictCox(m.cox, newdata = d[1:5,],
#'                    times = 1:5, type = "survival")
#' pred
#' 
#' pred.data <- predictCox(m.cox, newdata = d[1:4,],se=1L,
#'              times = 1:5, type = "survival", keep.newdata = TRUE)
#' pred.data
#'
#' m.cox <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X6,
#'                data = d, x = TRUE, y = TRUE)
#' pred.cox <- predictCox(m.cox, newdata = d[c(1:5,10,50),],
#'                        time = 1:5)
#' pred.cox
#' 
#' m.cox <- cph(Surv(time,event)~ strat(X1) + strat(X2) + X3 + X6,
#'              data = d, x = TRUE, y = TRUE)
#' pred.cox <- predictCox(m.cox, newdata = d[c(1:5,10,50),],
#'                        time = 1:5)
#' pred.cox
#' 
#' pred.dataci <- predictCox(m.cox, newdata = d[1:5,],
#'                        times = 1:5, keep.newdata = TRUE, se = TRUE)
#' pred.dataci
#'
#' @method print predictCox
#' @export
print.predictCox <- function(x,
                             digits = 3, ...){
    out <- as.data.table(x)
    print(out,digits=digits,...)
    invisible(out)
}



#----------------------------------------------------------------------
### print.predictCox.R ends here
