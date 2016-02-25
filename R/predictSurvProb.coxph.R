#' @title Predicting hazard or cumulative hazard
#' 
#' @aliases predictHazard predictHazard.coxph predictSurvProb predictSurvProb.coxph
#
#' @param object The fitted coxph model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param cause The event of interest
#' @param type should the hazard or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' Strange behavior of predict when strata
#' What to return if stratified Cox model with no covariate
#' 
#' @return a data table containing the predictions for each patient (in rows) and each time (in columns) and the strata (if any).
#' 
#' @examples 
#' 
#' library(data.table)
#' library(survival)
#' library(prodlim)
#' library(pec)    # needed to define the  predictSurvProb method
#' 
#' set.seed(10)
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- predictHazard.coxph(fit, newdata = d, times = 10)
#' res2 <- predictHazard.coxph(fit, newdata = d, times = d$time)
#' res3 <- predictSurvProb.coxph(fit, newdata = d, times = 10)
#' res4 <- predictSurvProb.coxph(fit, newdata = d, times = d$time)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- predictHazard.coxph(fitS, newdata = d, times = 10)
#' res2S <- predictHazard.coxph(fitS, newdata = d, times = d$time)
#' res3S <- predictSurvProb.coxph(fitS, newdata = d, times = 10)
#' res4S <- predictSurvProb.coxph(fitS, newdata = d, times = d$time)
#' 
#' @export
#' 

predictHazard.coxph <- function (object, newdata, times, cause = 1, type = "hazard",
                                 method.baseHaz = "dt") {
  
  require(data.table)
  match.arg(type, choices = c("hazard", "cumHazard"), several.ok = FALSE)
  
  times <- sort(times)
  Lambda0 <- baseHaz(object, cause = cause, method = method.baseHaz, centered = TRUE, lasttime =  times[length(times)],
                     addFirst = TRUE, addLast = TRUE) 
  
  strataspecials <- attr(object$terms,"specials")$strata 
  
  if(is.null(strataspecials)){
    Xb <- predict(object, newdata, type = "lp")
    
    ## compute hazard
    time.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times)
    resPred <- data.table( Lambda0[time.index, exp(Xb) %o% .SD[[1]], .SDcols = type] )
    setnames(resPred, paste0("t",times))
    
  }else{
    
    if("XXstrata" %in% names(newdata)){
      stop("predictSurvProb2.coxph: \'newdata\' must not contains a column named \"XXstrata\" \n")
    }
    if(!is.data.table(newdata)){
    newdata <- as.data.table(newdata)
    }
    
    ## linear predictor
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    Terms <- attr(object$terms, "term.labels")
    BaseVar <- Terms[ Terms %in% stratalabels == FALSE]
    if(length(BaseVar)>0){ # predict crashes if no additional variable
      Xb <-  rowSums(predict(object, newdata, type = "terms")) 
      # predict(object, newdata, type = "lp") has an output that I do not understand and that does not match "terms" in presence of strata
    }else{
      Xb <- rep(0, nrow(newdata))
    }
    
    ## define the strata
    stratalabels <- attr(object$terms,"term.labels")[strataspecials - 1]
    names.newdata<- names(newdata)
    strataVar <- names.newdata[which(paste0("strata(", names.newdata,")") %in% stratalabels)]
    nStrata <- length(strataVar)
    sapply(1:nStrata, function(x){
      newdata[, strataVar[x] := paste0(strataVar[x],"=",.SD[[1]]),.SDcols = strataVar[x], with = FALSE]
    })
    newdata[, XXstrata := interaction(.SD, drop = TRUE, sep = ", ") ,.SDcols = strataVar]
    
    ## compute hazard
    resPred <- NULL
    levelsStrata <- levels(newdata$XXstrata)
    for(iterS in levelsStrata){
      indexHaz <- Lambda0[,.I[strata == iterS]]
      indexNew <- newdata[,.I[XXstrata == iterS]]
      time.index <- prodlim::sindex(jump.times=Lambda0[indexHaz,time],eval.times=times)  
      resPred_tempo <- Lambda0[indexHaz[time.index], exp(Xb[indexNew]) %o% .SD[[1]], .SDcols = type]
      
      resPred <-  rbindlist(list(resPred, 
                                 data.table(resPred_tempo, strata = iterS, index = indexNew))
      )
    }
    setkey(resPred,index)
    resPred[,index := NULL]
    setnames(resPred, c(paste0("t",times),"strata"))
    
  }
  
  # export
  return(resPred)
}

predictSurvProb.coxph <- function(object, newdata, times, cause = 1,
                                  method.baseHaz = "dt") {
  
  res <- predictHazard.coxph(object, newdata = newdata, times = times, cause = cause, type = "cumHazard",
                             method.baseHaz = method.baseHaz)
  n.times <- length(times)
  
  ## convert cumulative hazard to survival
  res[, 1:n.times := exp(-.SD), .SDcols = 1:n.times]
 
  ## export
  return(res)
  
}

