#' @title Compute the baseline hazard
#
#' @param coxph the fitted coxph model
#' @param event the type of event corresponding to each observation [possible improvement: may be directly extracted from coxph ??]
#' @param time the time at which the event occured for each observation [possible improvement: may be directly extracted from coxph ??]
#' @param cause the event of interest
#' @param method the implementation to be used: "dt" or "cpp"
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' @return a data table containing the time, the strata (if any), the hazard and the cumulative hazard
#' 
#' @examples 
#' 
#' library(data.table)
#' library(rbenchmark)
#' library(prodlim)
#' library(rms)
#' library(testthat)
#' 
#' set.seed(10)
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- baseHaz(fit, d$status, d$time, method = "dt")
#' res2 <- baseHaz(fit, d$status, d$time, method = "cpp")
#' 
#' res3 <- baseHaz(fit, d$status, d$time, method = "dt", lasttime = 5)
#' res4 <- baseHaz(fit, d$status, d$time, method = "cpp", lasttime = 5)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- baseHaz(fitS, d$status, d$time, method = "dt")
#' res2S <- baseHaz(fitS, d$status, d$time, method = "cpp")
#' 
#' res3S <- baseHaz(fitS, d$status, d$time, method = "dt", lasttime = 5)
#' res4S <- baseHaz(fitS, d$status, d$time, method = "cpp", lasttime = 5)
#' 
#' @export
baseHaz <- function(object, cause = 1, method, center = TRUE, lasttime = Inf){
  
  # as in survival:::agsurv
  event <- object$y[, ncol(object$y)]
  time <- object$y[, ncol(object$y) - 1]
  
  ## strata
  strataspecials = attr(object$terms,"specials")$strata
  if (!is.null(strataspecials)){
    stratalabels = attr(object$terms,"term.labels")[strataspecials - 1]
    mod <- interaction(model.frame(object)[,stratalabels], drop = TRUE, sep = ", ")
  }
  
  ## method
  match.arg(method, choices = c("dt","cpp"), several.ok = FALSE)
  
  ## linear predictor
  lp <- object$linear.predictors
  if(center){
    lp <- lp + sum(object$means*coef(object))
  }
  
  ## main
  if (method == "cpp") {
    
    if (is.null(strataspecials)){
      resCpp <- BaseHazStrata_cpp(alltimes = time, status = event, Xb = lp, strata = NA,
                                  nPatients = object$n, nStrata = 1, lasttime = lasttime, cause = 1)
      
      resCpp$strata <- NULL
      
      if(is.na(resCpp$time[1])){ # failure or no event before lasttime
        resCpp <- matrix(nrow = 0,ncol = 3)
        colnames(resCpp) <- c("time", "hazard", "cumHazard")
      }
      
    }else{
      levelsStrata <- levels(mod)
      nStrata <- length(levelsStrata)
      
      resCpp <- BaseHazStrata_cpp(alltimes = time, status = event, Xb = lp, strata = as.numeric(mod) - 1,
                                    nPatients = object$n, nStrata = nStrata, lasttime = lasttime, cause = 1)
     
      if(is.na(resCpp$time[1])){ # failure or no event before lasttime
        resCpp <- matrix(nrow = 0,ncol = 4)
        colnames(resCpp) <- c("time", "strata", "hazard", "cumHazard")
      }else{
      resCpp$strata <- factor(resCpp$strata, levels = 0:(nStrata-1), labels = levelsStrata)
      }
    }
    
    return(as.data.table(resCpp))
    
  }else{ # method == "dt"
    
    #  R(ti) = {j ; Tj >= ti}
    #  sum_{j in R(ti)} sum_{k Tk = Tj} exp(beta^T Zk)
    require(data.table)
    
    dt.d <- data.table(time = time, event = event, lp = lp)
    if (is.null(strataspecials)){
      dt.d[,strata:=1]
      setorder(dt.d, time, -event)
    } else{
      dt.d[,strata:=mod]
      setorder(dt.d, strata,time, - event)
    }
    
    dt.d[, utime := cumsum(!duplicated(time)),by=strata]
    dt.d[, di := sum(event == cause), by = list(strata,utime)]
    dt.d[, elp := exp(lp)]
    dt.d[, Wti := sum(elp), by = list(strata,utime)]
    setkey(dt.d,strata,time)
    dt.d <- dt.d[unique(dt.d[,c("strata","time"),with=FALSE]),,mult="first"]
    ## dt.d <- dt.d[!duplicated(utime)]
    dt.d[, W := rev(cumsum(rev(Wti))),by=strata]
    
    if(!is.infinite(lasttime)){
      dt.d <- dt.d[time <= lasttime]
    }
    
      dt.d[, hazard := di/W]
      dt.d[, cumHazard := cumsum(hazard),by=strata]
    
      dt.d[,c("event","W","Wti","lp","elp","di","utime"):=NULL]
      
    if (is.null(strataspecials)) 
      dt.d[,strata:=NULL]
      
    return(dt.d)
  }
}
