#' @title Computing baseline hazard
#' 
#' @param object The fitted coxph or cph model
#' @param method The implementation to be used: "dt" or "cpp"
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
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- baseHazRR(fit, method = "dt")
#' res2 <- baseHazRR(fit, method = "cpp")
#' 
#' res3 <- baseHazRR(fit, method = "dt", maxtime = 5)
#' res4 <- baseHazRR(fit, method = "cpp", maxtime = 5)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- baseHazRR(fitS, method = "dt")
#' res2S <- baseHazRR(fitS, method = "cpp")
#' 
#' res3S <- baseHazRR(fitS, method = "dt", maxtime = 5)
#' res4S <- baseHazRR(fitS, method = "cpp", maxtime = 5)
#' @export
baseHazRR <- function(object,
                      method,
                      centered = TRUE,
                      maxtime = Inf){
    
    ## strata
    if ("cph"%in%class(object)){
        if(is.null(object$y)){
            stop("baseHazRR: argument \'y\' must be set to TRUE in cph \n")
        }
        strataspecials = attr(object$terms,"specials")$strat
    } else{
        strataspecials = attr(object$terms,"specials")$strata
    }
    test.Nstrata <- is.null(strataspecials)
    if (test.Nstrata == FALSE){
        stratavars = attr(object$terms,"term.labels")[strataspecials - 1]
        strataF <- interaction(model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
    }else{
        strataF <- factor("1")
    }
    
    ## method
    match.arg(method, choices = c("dt","cpp"), several.ok = FALSE)
    if(object$method == "exact"){
        stop("baseHazRR: exact correction for ties is not implemented \n")
    }
    if (is.na(maxtime)) maxtime <- Inf
    ## main
    time = object$y[, ncol(object$y) - 1]
    event = object$y[, ncol(object$y)]
    lp = if(centered == FALSE){object$linear.predictors + sum(object$means*coef(object))}else{object$linear.predictors}
    nPatients = length(time)
    method.ties = object$method
    method.baseHazRR = method
    maxtime = maxtime
    cause <- 1
    levelsStrata <- levels(strataF)
    nStrata <- length(levelsStrata)
    resCpp <- BaseHazStrata_cpp(alltimes = time,
                                status = event,
                                Xb = lp,
                                strata = as.numeric(strataF) - 1,
                                nPatients = nPatients,
                                nStrata = nStrata,
                                maxtime = maxtime,
                                cause = cause,
                                Efron = (method.ties == "efron"))

    if(length(resCpp$time) == 1 && is.na(resCpp$time[1])){ # failure or no event before maxtime
        warning("baseHazRR: failure \n")
        resCpp <- NULL
    }else{
        if (test.Nstrata) resCpp$strata <- factor(resCpp$strata, levels = 0:(nStrata-1), labels = levelsStrata)
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
