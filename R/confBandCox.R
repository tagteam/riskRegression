### confBandCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 21 2017 (17:40) 
## Version: 
## last-updated: apr 21 2017 (17:40) 
##           By: Brice Ozenne
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Compute confidence bands for a Cox model
#' @description Compute confidence bands for a Cox model
#' 
#' @param object
#' @param newdata
#' @param times
#' @param n.sim
#' @param alpah
#' 
#' @examples 
#' 
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' fit <- cox.aalen(Surv(time, event) ~ prop(X1) + prop(X2), data = d)
#' fit.coxph <- coxph(Surv(time, event) ~ X1 + X2, data = d, x = TRUE, y = TRUE)
#' newdata <- d[1:10]
#' 
#' set.seed(10)
#' res <- confBandCox(fit.coxph, newdata = newdata, n.sim = 500)
#' 
#' set.seed(10)
#' poutREF <- timereg:::predict.timereg(fit, newdata = newdata, times = fit$time.sim.resolution, resample.iid = 1)
#' 
#' # bug to be fixed: survival function should be right continuous
#' rbind(predictCox(fit.coxph, newdata = newdata[1], times = fit$time.sim.resolution-1e-8)$survival,
#'       predictCox(fit.coxph, newdata = newdata[1], times = fit$time.sim.resolution)$survival)
#' 
#' rbind(timereg:::predict.timereg(fit, newdata = newdata[1], times = fit$time.sim.resolution-1e-8, resample.iid = 0)$S0,
#'       timereg:::predict.timereg(fit, newdata = newdata[1], times = fit$time.sim.resolution, resample.iid = 0)$S0)
#' res
#' poutREF$unif.band


confBandCox <- function(object, newdata, times = NULL, n.sim, alpha = 0.05){
  # object <- fit.coxph
  # n.sim <- 500
  #
  #
  
  
  object.n <- CoxN(object)
  object.design <- CoxDesign(object)
  object.time <- c(0,unique(sort(object.design[object.design$status==1,"stop"])))-1e-5
  if(is.null(times)){
    times <- object.time
  }
  nTime <- length(times)
  
  new.n <- NROW(newdata)
  
  pred <- predictCox(object, newdata = newdata, times = times, se = TRUE, iid = TRUE, type = "survival")
  
  RR.se <- as.double(pred$survival.se/pred$survival)
  RR.delta <- sweep(pred$survival.iid, MARGIN = 1:2, FUN = "/", STATS = pred$survival)
  RR.delta <- matrix(RR.delta[], nrow = nTime*new.n, ncol = object.n)
  
  set.seed(10)
  RR.mpt <- .C("confBandBasePredict", delta = as.double(RR.delta), 
            nObs = as.integer(new.n), nt = as.integer(nTime), n = as.integer(object.n), 
            se = as.double(RR.se), mpt = double(n.sim * new.n), nSims = as.integer(n.sim), 
            PACKAGE = "timereg")$mpt
  RR.mpt <- matrix(RR.mpt, n.sim, new.n, byrow = TRUE)
  RR.uband <- apply(RR.mpt, 2, percen, per = 1 - alpha)
  
  return(RR.uband)
}


#----------------------------------------------------------------------
### confBandCox.R ends here
