#' @title Extrating predicting risks with confidence intervals from regression models 
#' 
#' @description Extract event probabilities with confidence intervals from fitted regression models and machine learning objects. 
#' 
#' The function predictRiskCI is a generic function, meaning that it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' See \code{\link{predictRiskCI}}.
#' 
#' @aliases predictRiskCI predictRiskCI.CauseSpecificCox
#' predictRiskCI.coxph predictRiskCI.cph
#' 
#' @usage
#' \method{predictRiskCI}{cph}(object,newdata,conf.level=0.95,times,colnames=TRUE,...)
#' \method{predictRiskCI}{coxph}(object,newdata,conf.level=0.95,times,colnames=TRUE,...)
#' \method{predictRiskCI}{CauseSpecificCox}(object,newdata,conf.level=0.95,times,cause,colnames=TRUE,...)
#' @inheritParams predictRisk
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param conf.level confidence level of the interval.
#' @param colnames Logical. If \code{TRUE} set the column names to he evaluation times in the output.
#' @return A list containing the punctual estimate, the lower bound and the upper bound of the confidence interval.
#' Each element of the list is organised in the same way as the output of \code{predictRisk}.
#' @author  Brice Ozenne broz@@sund.ku.dk
#' @seealso See \code{\link{predictRisk}}.
#' @keywords survival
#' @examples
#' set.seed(100)
#' 
#' #### survival ####
#' d <- sampleData(200,outcome="survival")
#' # fit a Cox model
#' library(rms)
#' cphmodel <- cph(Surv(time,event)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
#' library(survival)
#' coxphmodel <- coxph(Surv(time,event)~X1+X2,data=d,x=TRUE,y=TRUE)
#'
#' # prediction with confidence interval
#' predictRiskCI(cphmodel, newdata = d[1:5,], times = 1:3)
#' predictRiskCI(coxphmodel, newdata = d[1:5,], times = 1:3)
#'
#' #### competing risks ####
#' d <- sampleData(200,outcome="competing.risks")
#' CSCmodel <- CSC(Hist(time,event)~X1+X2, data = d)
#'
#' predictRiskCI(CSCmodel, newdata = d[1:5,], times = 1:3, cause = 1)
#' @export 
predictRiskCI <- function(object, newdata, conf.level,...){
  UseMethod("predictRiskCI",object)
}

##' @export
predictRiskCI.coxph <- function(object,
                                newdata,                                
                                conf.level=0.95,
                                times,
                                colnames=TRUE,
                                ...){

   res <- predictCox(object=object,
                      newdata = newdata,
                      times = times,
                      se = TRUE,
                      iid = FALSE,
                      keep.times=colnames,
                      keep.lastEventTime=FALSE,
                      type="survival")

    p <- list(estimate = 1-res$survival)
    p$lower <- p$estimate - 1.96 * res$survival.se
    p$upper <- p$estimate + 1.96 * res$survival.se

    if(colnames){
        colnames(p$estimate) <- res$times
        colnames(p$lower) <- res$times
        colnames(p$upper) <- res$times
    }
        
    if (NROW(p$estimate) != NROW(newdata) || NCOL(p$estimate) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p$estimate)," x ",NCOL(p$estimate),"\n\n",sep=""))
    }

    if(any(p$lower<0) || any(p$upper>1)){
        warning("Some of the bounds of the confidence intervals are not in [0;1]. \n Wald-type confidence intervals may not be appropriate for all predictions \n")
    }
    return(p)
}

##' @export
predictRiskCI.cph <- predictRiskCI.coxph

##' @export
predictRiskCI.CauseSpecificCox <- function(object,
                                           newdata,
                                           conf.level = 0.95,
                                           times,
                                           cause,
                                           colnames = TRUE,
                                           ...) { 
    res <- predict(object = object,
                 newdata = newdata,
                 times = times,
                 cause = cause,
                 keep.strata = FALSE,
                 keep.times = colnames,
                 se = TRUE,
                 iid = FALSE)
        
    p <- list(estimate = res$absRisk,
              lower = res$absRisk - 1.96 * res$absRisk.se,
              upper = res$absRisk + 1.96 * res$absRisk.se
              )

    if(colnames){
        colnames(p$estimate) <- res$times
        colnames(p$lower) <- res$times
        colnames(p$upper) <- res$times
    }
    
    if (NROW(p$estimate) != NROW(newdata) || NCOL(p$estimate) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p$estimate)," x ",NCOL(p$estimate),"\n\n",sep=""))
    }

    if(any(p$lower<0)){
        warning("Some of the bounds of the confidence intervals are below 0. \n Wald-type confidence intervals may not be appropriate for all predictions \n")
    }
    return(p)

}



