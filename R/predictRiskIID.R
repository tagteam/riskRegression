#' @title Extrating the influence function from regression models 
#' 
#' @description Extract the influence function from fitted regression models and machine learning objects. 
#' 
#' The function predictRiskCI is a generic function, meaning that it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' See \code{\link{predictRiskIID}}.
#' 
#' @aliases predictRiskIID predictRiskIID.CauseSpecificCox
#' predictRiskIID.coxph predictRiskIID.cph
#' 
#' @usage
#' \method{predictRiskIID}{cph}(object,newdata,times,...)
#' \method{predictRiskIID}{coxph}(object,newdata,times,...)
#' \method{predictRiskIID}{CauseSpecificCox}(object,newdata,times,cause,...)
#' @inheritParams predictRisk
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
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
#' predictRiskIID(cphmodel, newdata = d[1:5,], times = 1:3)
#' predictRiskIID(coxphmodel, newdata = d[1:5,], times = 1:3)
#'
#' #### competing risks ####
#' d <- sampleData(200,outcome="competing.risks")
#' CSCmodel <- CSC(Hist(time,event)~X1+X2, data = d)
#'
#' predictRiskIID(CSCmodel, newdata = d[1:5,], times = 1:3, cause = 1)
#' 
#' @export 
predictRiskIID <- function(object, newdata,...){
  UseMethod("predictRiskIID",object)
}

##' @export
predictRiskIID.coxph <- function(object,
                                 newdata,
                                 times,    
                                 ...){

    iid <- predictCox(object=object,
                      newdata = newdata,
                      times = times,
                      se = FALSE,
                      iid = TRUE,
                      keep.times=FALSE,
                      keep.lastEventTime=FALSE,
                      type="survival")$survival.iid

    return(iid)
}

##' @export
predictRiskIID.cph <- predictRiskIID.coxph

##' @export
predictRiskIID.CauseSpecificCox <- function (object,
                                             newdata,
                                             times,
                                             cause,
                                             ...){
              
    iid <- predict(object = object,
                   newdata = newdata,
                   times = times,
                   cause = cause,
                   keep.strata = FALSE,
                   keep.times = FALSE,
                   se = FALSE,
                   iid = TRUE)$absRisk.iid
       
    return(iid)

}



