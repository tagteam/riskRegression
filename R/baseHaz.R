#' @export
baseHaz <- function (x, ...) {
  UseMethod("baseHaz", x)
}

.CST_EPSILON <- 1e-10 # time gap after the last event for the survival to be set to NA

#' @title Computing baseline hazard
#' 
#' @aliases baseHaz baseHaz.coxph
#' 
#' @param object The fitted coxph model
#' @param method The implementation to be used: "dt" or "cpp"
#' @param centered If TRUE remove the centering factor used by coxph in the linear predictor
#' @param lasttime Baseline hazard will be computed for each event before lasttime
#' @param addFirst Should the hazard at time 0 be also returned
#' @param addLast Should the hazard after the last event be also returned
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
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- baseHaz(fit, method = "dt")
#' res2 <- baseHaz(fit, method = "cpp")
#' 
#' res3 <- baseHaz(fit, method = "dt", lasttime = 5)
#' res4 <- baseHaz(fit, method = "cpp", lasttime = 5)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- baseHaz(fitS, method = "dt")
#' res2S <- baseHaz(fitS, method = "cpp")
#' 
#' res3S <- baseHaz(fitS, method = "dt", lasttime = 5)
#' res4S <- baseHaz(fitS, method = "cpp", lasttime = 5)
#' 
#' @export
baseHaz.coxph <- function(object, method, centered = TRUE, lasttime = Inf, addFirst = FALSE, addLast = FALSE){
  
  # as in survival:::agsurv
  cause <- 1
  event <- object$y[, ncol(object$y)]
  time <- object$y[, ncol(object$y) - 1]
  
  ## strata
  strataspecials = attr(object$terms,"specials")$strata
  if (!is.null(strataspecials)){
    stratalabels = attr(object$terms,"term.labels")[strataspecials - 1]
    mod <- interaction(model.frame(object)[,stratalabels], drop = TRUE, sep = ", ", lex.order = TRUE) # [COULD BE OPTIMIZED in C]
  }else{
    mod <- factor("1")
  }
  
  ## method
  match.arg(method, choices = c("dt","cpp"), several.ok = FALSE)
  if(object$method == "exact"){
    stop("baseHaz.coxph: exact correction for ties is not implemented \n")
  }
  
  ## linear predictor
  lp <- object$linear.predictors
  if(centered == FALSE){
    lp <- lp + sum(object$means*coef(object))
  }
  
  ## main
  if (method == "cpp") {
    
    if (is.null(strataspecials)){
      resCpp <- BaseHazStrata_cpp(alltimes = time, status = event, Xb = lp, strata = NA,
                                  nPatients = object$n, nStrata = 1, lasttime = lasttime, cause = cause,
                                  Efron = (object$method == "efron"), addFirst = addFirst, addLast = addLast)
      
      resCpp$strata <- NULL
      
      if(length(resCpp$time) == 1 && is.na(resCpp$time[1])){ # failure
        warning("baseHaz.coxph: failure \n")
        resCpp <- matrix(nrow = 0,ncol = 3)
        colnames(resCpp) <- c("time", "hazard", "cumHazard")
      }
      
    }else{
      levelsStrata <- levels(mod)
      nStrata <- length(levelsStrata)
      
      resCpp <- BaseHazStrata_cpp(alltimes = time, status = event, Xb = lp, strata = as.numeric(mod) - 1,
                                  nPatients = object$n, nStrata = nStrata, lasttime = lasttime, cause = cause,
                                  Efron = (object$method == "efron"), addFirst = addFirst, addLast = addLast)
      
      if(length(resCpp$time) == 1 && is.na(resCpp$time[1])){ # failure or no event before lasttime
        warning("baseHaz.coxph: failure \n")
        resCpp <- matrix(nrow = 0,ncol = 4)
        colnames(resCpp) <- c("time", "strata", "hazard", "cumHazard")
      }else{
        resCpp$strata <- factor(resCpp$strata, levels = 0:(nStrata-1), labels = levelsStrata)
      }
    }
    
    return(data.table::as.data.table(resCpp))
    
  }else{ # method == "dt"
    
    #  R(ti) = {j ; Tj >= ti}
    #  sum_{j in R(ti)} sum_{k Tk = Tj} exp(beta^T Zk)
    
    dt.d <- data.table::data.table(time = time, event = event, lp = lp)
    
    ## strata
    if (is.null(strataspecials)){
      dt.d[,strata:="1"]
      setorder(dt.d, time, -event)
    } else{
      dt.d[,strata:=mod]
      setorder(dt.d, strata,time, - event)
    }
    
    ## individual level
    dt.d[, utime := cumsum(!duplicated(time)),by=strata]
    dt.d[, di := sum(event == cause), by = list(strata,utime)]
    dt.d[, elp := exp(lp)]
    dt.d[, Wti := sum(elp), by = list(strata,utime)]
    if(object$method == "efron"){
      dt.d[, Wti_event := sum(elp*(event == cause)), by = list(strata,utime)]
    }
    data.table::setkey(dt.d,strata,time)
    
    ## group level
    dt.d <- dt.d[unique(dt.d[,c("strata","time"),with=FALSE]),,mult="first"]
    dt.d[, W := rev(cumsum(rev(Wti))),by=strata]
    
    if (addLast){
      maxTimeEvent <- dt.d[, time[.N], by = strata][[2]]
    }
    if(!is.infinite(lasttime)){
      dt.d <- dt.d[time <= lasttime]
    }
    
    ## compute hazard
    if(object$method == "efron"){
      dt.d[, hazard := di * baseHazEfron_survival_cpp(ntimes = .N, ndead = di, risk = W, riskDead = Wti_event)]
      dt.d[,Wti_event := NULL]
    }else if(object$method == "breslow"){
      dt.d[, hazard := di/W]
    }
    
    dt.d[, cumHazard := cumsum(hazard), by = strata]
    
    ## remove useless variables
    dt.d[,c("event","W","Wti","lp","elp","di","utime") := NULL]
    
    ## add before or after first or last event
    if (addFirst){ 
      dt.d <- dt.d[, .(time = c(0,time), hazard = c(0,hazard), cumHazard = c(0,cumHazard))  , by = strata]
      if( length(unique(dt.d$strata)) != length(levels(mod)) ){ # if some strata has been removed due to lasttime
        dt.d <- data.table::rbindlist(list(dt.d,
                                           data.table::data.table(strata = setdiff(levels(mod),unique(dt.d$strata)), time = 0, hazard = 0, cumHazard = 0)))
        data.table::setkey(dt.d,strata,time)
      }
      data.table::setcolorder(dt.d, c("time", "strata", "hazard", "cumHazard"))
    }
    
    if (addLast){
      dt.d <- dt.d[, .(time = c(time, maxTimeEvent[.GRP] + .CST_EPSILON), hazard = c(hazard,NA), cumHazard = c(cumHazard,NA))  , by = strata]
      
      if( length(unique(dt.d$strata)) != length(levels(mod)) ){# if some strata has been removed due to lasttime
        index.diff <- which(levels(mod) %in% unique(dt.d$strata) == FALSE)
        dt.d <- data.table::rbindlist(list(dt.d,
                                           data.table::data.table(strata = levels(mod)[index.diff], time = maxTimeEvent[index.diff] + .CST_EPSILON, hazard = NA, cumHazard = NA)))
        data.table::setkey(dt.d,strata,time)
      }
      data.table::setcolorder(dt.d, c("time", "strata", "hazard", "cumHazard"))
    }
    
    ## remove strata
    if (is.null(strataspecials)){ 
      dt.d[,strata:=NULL]
    }
    
    ## export
    return(dt.d)
  }
}
