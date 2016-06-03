#' @title Computing baseline hazard
#'
#' Fast routine to get baseline hazards and subject specific hazards as well as survival probabilities from a coxph or cph object
#' @param object The fitted coxph or cph model
#' @param centered If TRUE remove the centering factor used by coxph in the linear predictor
#' @param maxtime Baseline hazard will be computed for each event before maxtime
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' @return a data table containing the time, the strata (if any), the hazard and the cumulative hazard.
#' 
#' @examples 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#' 
#' baseHazRR(fit)
#' cbind(survival::basehaz(fit),baseHazRR(fit))
#' baseHazRR(fit, maxtime = 5)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- baseHazRR(fitS)
#' res2S <- baseHazRR(fitS, maxtime = 5)
#' @export
baseHazRR <- function(object,
                      centered = TRUE,
                      maxtime = Inf){
    
  ## extract elements from objects
  if ("cph" %in% class(object)){
    
    nPatients <- sum(object$n)
    if(is.null(object$y)){
      stop("baseHazRR: argument \'y\' must be set to TRUE in cph \n")
    }
    strataspecials <- attr(object$terms,"specials")$strat
    no.strata <- is.null(strataspecials)
    if(!no.strata){
      strataF <- object$Strata  
    }else{
      strataF <- factor("1")
    }
    
    
  } else if ("coxph" %in% class(object)){
    
      nPatients <- object$n
      strataspecials <- attr(object$terms,"specials")$strata
      no.strata <- is.null(strataspecials)
      if(!no.strata){
          stratavars <- attr(object$terms,"term.labels")[strataspecials - 1]
          strataF <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ".", lex.order = TRUE) 
      }else{
          strataF <- factor("1")
      }
    
  } else {
      stop("baseHazRR: only implemented for \"coxph\" and \"cph\" objects \n")
  }
  
  
  ## method
  if(object$method == "exact"){
    stop("baseHazRR: exact correction for ties is not implemented \n")
  }
  if (is.na(maxtime)) maxtime <- Inf
  
    ## main
    levelsStrata <- levels(strataF)
    nStrata <- length(levelsStrata)
  
  resCpp <- BaseHazStrata_cpp(alltimes = object$y[, ncol(object$y) - 1],
                              status = object$y[, ncol(object$y)],
                              Xb = if(centered == FALSE){object$linear.predictors + sum(object$means*coef(object))}else{object$linear.predictors},
                              strata = as.numeric(strataF) - 1,
                              nPatients = nPatients,
                              nStrata = nStrata,
                              maxtime = maxtime,
                              cause = 1,
                              Efron = (object$method == "efron"))
  
  if (no.strata == FALSE){ ## rename the strata value with the correct levels
    resCpp$strata <- factor(resCpp$strata, levels = 0:(nStrata-1), labels = levelsStrata)
  }
  
  return(data.table::as.data.table(resCpp))
}

# else{ # method.baseHazRR == "dt"

#  R(ti) = {j ; Tj >= ti}
#  sum_{j in R(ti)} sum_{k Tk = Tj} exp(beta^T Zk)

## dt.d <- data.table::data.table(time = time, event = event, lp = lp)

## strata
## if (test.Nstrata){
## dt.d[,strata:="1"]
## setorder(dt.d, time, -event)
## levelsStrata <- "1"
## } else{
## dt.d[,strata := strataF]
## setorder(dt.d, strata,time, - event)
## levelsStrata <- levels(strataF)
## }

## individual level
## dt.d[, utime := cumsum(!duplicated(time)),by=strata]
## dt.d[, di := sum(event == cause), by = list(strata,utime)]
## dt.d[, elp := exp(lp)]
## dt.d[, Wti := sum(elp), by = list(strata,utime)]
## if(method.ties == "efron"){
## dt.d[, Wti_event := sum(elp*(event == cause)), by = list(strata,utime)]
## }
## data.table::setkey(dt.d,strata,time)

## group level
## dt.d <- dt.d[unique(dt.d[,c("strata","time"),with=FALSE]),,mult="first"]
## dt.d[, W := rev(cumsum(rev(Wti))),by=strata]

## if(!is.infinite(maxtime)){
## dt.d <- dt.d[time <= maxtime]
## }

## compute hazard
## if(method.ties == "efron"){
## dt.d[, hazard := di * baseHazEfron_survival_cpp(ntimes = .N, ndead = di, risk = W, riskDead = Wti_event)]
## dt.d[,Wti_event := NULL]
## }else if(method.ties == "breslow"){
## dt.d[, hazard := di/W]
## }

## dt.d[, cumHazard := cumsum(hazard), by = strata]

## remove obsolete variables
## dt.d[,c("event","W","Wti","lp","elp","di","utime") := NULL]

## remove strata
## if (test.Nstrata){ 
## dt.d[,strata:=NULL]
## }

## export
## return(dt.d)
## }

## }
