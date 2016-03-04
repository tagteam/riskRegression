baseHaz <- function (x, ...) {
  UseMethod("baseHaz", x)
}


#' @title Computing baseline hazard
#
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
#' When considering strata, return the result in a different order compared to basehaz
#' interaction not optimal in term of computation time
#' 
#' @return a data table containing the time, the strata (if any), the hazard and the cumulative hazard.
#' 
#' @examples 
#' 
#' library(data.table)
#' library(survival)
#' library(prodlim)
#' 
#' set.seed(10)
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- baseHaz(fit, method = "dt", centered = FALSE)
#' res2 <- baseHaz(fit, method = "cpp", centered = FALSE)
#' 
#' res3 <- baseHaz(fit, method = "dt", centered = FALSE, lasttime = 5)
#' res4 <- baseHaz(fit, method = "cpp", centered = FALSE, lasttime = 5)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- baseHaz(fitS, method = "dt", centered = FALSE)
#' res2S <- baseHaz(fitS, method = "cpp", centered = FALSE)
#' 
#' res3S <- baseHaz(fitS, method = "dt", centered = FALSE, lasttime = 5)
#' res4S <- baseHaz(fitS, method = "cpp", centered = FALSE, lasttime = 5)
#' 
#' @export
#' 
baseHaz.coxph <- function(object, method, centered = TRUE, lasttime = Inf, addFirst = FALSE, addLast = FALSE){

  # as in survival:::agsurv
  cause <- 1
  event <- object$y[, ncol(object$y)]
  time <- object$y[, ncol(object$y) - 1]
  
  ## strata
  strataspecials = attr(object$terms,"specials")$strata
  if (!is.null(strataspecials)){
    stratalabels = attr(object$terms,"term.labels")[strataspecials - 1]
    mod <- interaction(model.frame(object)[,stratalabels], drop = TRUE, sep = ", ") # [COULD BE OPTIMIZED in C]
  }
  
  ## covariables
  if(object$method == "efron"){
    covlabels <- attr(object$terms,"term.labels")
    if (!is.null(strataspecials)){
      covlabels <- setdiff(covlabels, stratalabels)
      }
  }
  
  ## method
  match.arg(method, choices = c("dt","cpp"), several.ok = FALSE)
  if(object$method == "exact"){
    stop("baseHaz.coxph: exact correction for ties is not implemented \n")
  }else if(object$method == "efron" && method == "cpp"){
    stop("baseHaz.coxph: efron correction for ties is not (yet) implemented for method = \"cpp\" \n",
         "use method = \"dt\" instead \n")
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
                                  Efron = Efron, addFirst = addFirst, addLast = addLast)
      
      resCpp$strata <- NULL
      
      if(is.na(resCpp$time[1])){ # failure or no event before lasttime
        resCpp <- matrix(nrow = 0,ncol = 3)
        colnames(resCpp) <- c("time", "hazard", "cumHazard")
      }
      
    }else{
      levelsStrata <- levels(mod)
      nStrata <- length(levelsStrata)
      
      resCpp <- BaseHazStrata_cpp(alltimes = time, status = event, Xb = lp, strata = as.numeric(mod) - 1,
                                  nPatients = object$n, nStrata = nStrata, lasttime = lasttime, cause = cause,
                                  Efron = Efron, addFirst = addFirst, addLast = addLast)
     
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
    if(object$method == "efron"){
      dt.d[, covlabels := sweep(model.frame(coxph.E)[,covlabels], MARGIN = 2, FUN = "-", STATS = object$means), 
           with = FALSE]
    }
    ## strata
    if (is.null(strataspecials)){
      dt.d[,strata:=1]
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
      dt.d[, c(paste0("temp.",covlabels)) := lapply(covlabels, function(x){sum(.SD[[1]] * .SD[[x]])}), 
           .SDcols = c("elp",covlabels), 
           by = list(strata,utime)]
      dt.d[, c(paste0("xsum2.",covlabels)) := lapply(covlabels, function(x){sum(.SD[[1]] * .SD[[2]] * .SD[[x]])}), 
           .SDcols = c("elp","event",covlabels), 
           by = list(strata,utime)]
    }
    
    setkey(dt.d,strata,time)
     
    ## group level
    dt.d <- dt.d[unique(dt.d[,c("strata","time"),with=FALSE]),,mult="first"]
    dt.d[, W := rev(cumsum(rev(Wti))),by=strata]
    
    ## compute hazard
    if(object$method == "efron"){
      
      dt.d[, c(paste0("xsum.",covlabels)) := lapply(.SD, function(x){rev(cumsum(rev(x)))}),
          .SDcols = c(paste0("temp.",covlabels)), by=strata][,-1, with = FALSE]
      
      n.dt <- nrow(dt.d)
      n.cov <- length(covlabels)
      
#       dt.d[, hazard := di * baseHazEfron_survival_cpp(n= n.dt, nvar = n.cov, dd = .SD$di, x1 = .SD$W, x2 = .SD$Wti,
#                                                       xsum = unlist(.SD[,c(paste0("xsum.",covlabels)), with = FALSE]), 
#                                                       xsum2 = unlist(.SD[,c(paste0("xsum2.",covlabels)), with = FALSE])
#                                                       )$sum1, 
#            .SDcols = c("di","W","Wti", c(paste0("xsum.",covlabels), c(paste0("xsum2.",covlabels)))),
#            by=strata]
            dt.d[, hazard := di * baseHazEfron_survival_cpp(ntimes = n.dt, nvar = n.cov, ndead = .SD$di, risk = .SD$W, riskDead = .SD$Wti), 
                 .SDcols = c("di","W","Wti", c(paste0("xsum.",covlabels), c(paste0("xsum2.",covlabels)))),
                 by=strata]
      
       dt.d[, cumHazard := cumsum(hazard), by=strata]
      
      if(!is.infinite(lasttime)){
        dt.d <- dt.d[time <= lasttime]
      }
      
      dt.d[,c(paste0("temp.",covlabels),paste0("xsum.",covlabels),paste0("xsum2.",covlabels)):=NULL]
      
    }else if(object$method == "breslow"){
      
      if(!is.infinite(lasttime)){
        dt.d <- dt.d[time <= lasttime]
      }
      
      dt.d[, hazard := di/W]
      dt.d[, cumHazard := cumsum(hazard),by=strata]
      
    }
    
    ## remove useless variables
    dt.d[,c("event","W","Wti","lp","elp","di","utime"):=NULL]
      
    ## add before or after first or last event
    if (addFirst){ 
      dt.d <- dt.d[, .(time = c(0,time), hazard = c(0,hazard), cumHazard = c(0,cumHazard))  , by = strata]
    }
    if (addLast){
      dt.d <- dt.d[, .(time = c(time, time[.N] + 1e-10), hazard = c(hazard,NA), cumHazard = c(cumHazard,NA))  , by = strata]
      setcolorder(dt.d, c("time", "strata", "hazard", "cumHazard"))
    }
    ## remove strata
    if (is.null(strataspecials)){ 
      dt.d[,strata:=NULL]
    }

    ## export
    return(dt.d)
  }
}
