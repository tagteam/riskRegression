#' @title Compute the baseline hazard
#
#' @param coxph the fitted coxph model
#' @param event the type of event corresponding to each observation [possible improvement: may be directly extracted from coxph ??]
#' @param time the time at which the event occured for each observation [possible improvement: may be directly extracted from coxph ??]
#' @param cause the event of interest
#' @param method the implementation to be used: "dt" or "cpp"
#' 
#' @details 
#' "cpp" seems to be faster up to 10 000 observations. For very large number of observations "dt" seems faster [I do not understant why]
#' Because of the way of computing W, "cpp" is sensible to rounding error for instance when one coefficient is large compared to the others (then W -= exp(Xb) may be problematic)
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
#' #### check validity
#' res0 <- basehaz(fit, centered = FALSE)
#' res1 <- basehazFast(fit, d$status, d$time, method = "dt")
#' res2 <- basehazFast(fit, d$status, d$time, method = "cpp")
#' test_that("basehazFast",{
#'   expect_equal(res0, as.data.frame(res1))    
#' })
#' 
#' test_that("basehazFast",{
#'   expect_equal(res0, as.data.frame(res2))    
#' })
#' 
#' #### check timing
#' ben <- benchmark(replications=100,
#' basehaz(fit, centered = FALSE),
#' basehazFast(fit, d$status, d$time, method = "dt"),
#' basehazFast(fit, d$status, d$time, method = "cpp"),
#' columns=c('test', 'replications', 'elapsed',"user.self","sys.self")
#' )
#' ben
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' basehazFast(fitS, d$status, d$time, method = "dt")
#' 
#' @export
basehazFast <- function(object, event, time, cause = 1, method){

    ## strata
    strataspecials = attr(object$terms,"specials")$strata
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    strata <- object$xlevels[[stratalabels]] 
    mod <- model.frame(object) 
    Terms <- attr(object$terms, "term.labels")
    ## method
    match.arg(method, choices = c("dt","cpp"), several.ok = FALSE)
    require(data.table)
  
    dt.d <- data.table(time = time, event = event,
                       lp = object$linear.predictors + sum(object$means*coef(object)))
    if (is.null(strataspecials)){
        dt.d[,strata:=1]
        setorder(dt.d, time, -event)
    } else{
          dt.d[,strata:=mod[,stratalabels,drop=TRUE]]
          setorder(dt.d, strata,time, - event)
      }
    if(method == "cpp"){
        Lambda0_cpp <- BaselineHazardFast_cpp(status = dt.d$event, 
                                              pattimes = dt.d$time, 
                                              Xb = dt.d$lp,
                                              n_patients = object$n,
                                              eventtimes = unique(dt.d$time),
                                              cause = cause)
    
    return(data.table(hazard = Lambda0_cpp$hazard,
                      time = Lambda0_cpp$time))
  }else{ # method == "dt"
    
    #  R(ti) = {j ; Tj >= ti}
    #  sum_{j in R(ti)} sum_{k Tk = Tj} exp(beta^T Zk)
       
       dt.d[, utime := cumsum(!duplicated(time)),by=strata]
       dt.d[, di := sum(event == cause), by = list(strata,utime)]
       dt.d[, elp := exp(lp)]
       dt.d[, Wti := sum(elp), by = list(strata,utime)]
       setkey(dt.d,strata,time)
       dt.d <- dt.d[unique(dt.d[,c("strata","time"),with=FALSE]),,mult="first"]
       ## dt.d <- dt.d[!duplicated(utime)]
       dt.d[, W := rev(cumsum(rev(Wti))),by=strata]
       dt.d[, h := di/W]
       dt.d[, H := cumsum(h),by=strata]
       dt.d[,c("event","W","Wti","lp","elp","di","utime"):=NULL]
       if (is.null(strataspecials)) 
           dt.d[,strata:=NULL]
       return(dt.d)
       ## [,.(hazard = h,cumulative.hazard = H, time = time,strata=strata)])
   }
}
