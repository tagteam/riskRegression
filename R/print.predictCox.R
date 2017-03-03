### print.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 15 2017 (17:36) 
## Version: 
## last-updated: Mar  3 2017 (09:22) 
##           By: Thomas Alexander Gerds
##     Update #: 138
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
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param conf.level confidence level of the interval.
#' @param digit integer indicating the number of decimal places.
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
#' pred.data <- predictCox(m.cox, newdata = d[1:4,],
#'              times = 1:5, type = "survival", keep.newdata = TRUE)
#' pred.data
#' print(pred.data, reduce.data = TRUE)
#'
#' pred.ci <- predictCox(m.cox, newdata = d[1:5,],
#'                       times = 1:5,se = TRUE)
#' print(pred.ci, ci = TRUE)
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
#' print(pred.dataci, ci = TRUE)
#'
#' @method print predictCox
#' @export
print.predictCox <- function(x,
                             digits = 3, ...){
    if (is.null(x$newdata)){
        print.listof(x)
    } else{
        nd=x$newdata
        data.table::setDT(nd)
        out <- data.table::rbindlist(lapply(1:length(x$times),function(tt){
            ndtt=copy(nd)
            nd[,times:=x$times[[tt]]]
            if (!is.null(x$strata))
                nd[,strata:=x$strata]
            for (name in x$type){
                tyc <- cbind(x[[name]][,tt])
                colnames(tyc) <- name
                if (x$se==1L){
                    tyc <- cbind(tyc,x[[paste0(name,".se")]][,tt],x[[paste0(name,".lower")]][,tt],x[[paste0(name,".upper")]][,tt])
                    colnames(tyc) <- paste0(name,c("",".se",".lower",".upper"))
                }
                ## setDT(tyc)
                nd <- cbind(nd,tyc)
            }
            nd   
        }))
        print(out,digits=digits,...)
        invisible(out)
    }
}


#----------------------------------------------------------------------
### print.predictCox.R ends here
