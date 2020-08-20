### influence2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 19 2017 (16:01) 
## Version: 
## last-updated: aug 19 2020 (09:17) 
##           By: Brice Ozenne
##     Update #: 184
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Influence test (documentation)
#' @title Influence test [Experimental!!]
#' @description Compare two estimates using their influence function
#' @name influenceTest
#' 
#' @param object either a list of models or an object of class predictCox or predictCSC.
#' @param object2 same as predict1 but for another model.
#' @param type [character]the type of predicted value. 
#' @param newdata [data.frame or data.table] Contain
#'     the values of the predictor variables defining subject specific
#'     predictions.
#' @param times [numeric vector] Time points at which to return
#' the estimated absolute risk.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the covariates
#' used to make the prediction in the output. 
#' @param keep.strata [logical] If \code{TRUE} add the value of the strata
#' used to make the prediction in the output. 
#' @param cause [integer/character] Identifies the cause of interest among the competing
#'     events.
#' @param band [logical] If \code{TRUE} add the influence function to the output
#' such that \code{confint} will be able to compute the confidence bands. 
#' @param ... additional arguments to be passed to lower level functions.
#' 
#' @export
`influenceTest` <-
  function(object,...) UseMethod("influenceTest")

## * Influence test (code)
#' @examples
#' library(lava)
#' library(survival)
#' library(prodlim)
#' library(data.table)
#' n <- 100
#'
#' #### Under H1
#' set.seed(1)
#' newdata <- data.frame(X1=0:1)
#'
#' ## simulate non proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' regression(m) <- s ~ exp(-2*X1)
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- coxWeibull.lvm(scale=1,shape=~s)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#'
#' ## fit cox models
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, 
#'                data = d, y = TRUE, x = TRUE)
#'
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
#'                      data = d, y = TRUE, x = TRUE)
#'
#' ## compare models
#' # one time point
#' outIF <- influenceTest(list(m.cox, mStrata.cox), 
#'               type = "survival", newdata = newdata, times = 0.5)
#' confint(outIF)
#'                                  
#' # several timepoints
#' outIF <- influenceTest(list(m.cox, mStrata.cox), 
#'               type = "survival", newdata = newdata, times = c(0.5,1,1.5))
#' confint(outIF)
#'
#' #### Under H0 (Cox) ####
#' set.seed(1)
#' ## simulate proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm()
#' distribution(m,~y) <- coxWeibull.lvm()
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#' 
#' ## fit cox models
#' Utime <- sort(unique(d$eventtime))
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, 
#'                data = d, y = TRUE, x = TRUE)
#'
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
#'                      data = d, y = TRUE, x = TRUE)
#'
#' p.cox <- predictCox(m.cox, newdata = newdata, time = Utime, type = "survival")
#' p.coxStrata <- predictCox(mStrata.cox, newdata = newdata, time = Utime, type = "survival")
#'
#' ## display
#' library(ggplot2)
#' autoplot(p.cox)
#' autoplot(p.coxStrata)
#'  
#' ## compare models
#' outIF <- influenceTest(list(m.cox, mStrata.cox), 
#'                        type = "survival", newdata = newdata, times = Utime[1:6])
#' confint(outIF)
#' 
#' #### Under H0 (CSC) ####
#' set.seed(1)
#' ff <- ~ f(X1,2) + f(X2,-0.033)
#' ff <- update(ff, ~ .+ f(X3,0) + f(X4,0) + f(X5,0))
#' ff <- update(ff, ~ .+ f(X6,0) + f(X7,0) + f(X8,0) + f(X9,0))
#' d <- sampleData(n, outcome = "competing.risk", formula = ff)
#' d[,X1:=as.numeric(as.character(X1))]
#' d[,X2:=as.numeric(as.character(X2))]
#' d[,X3:=as.numeric(as.character(X3))]
#' d[,X4:=as.numeric(as.character(X4))]
#' d[,X5:=as.numeric(as.character(X5))]
#' setkey(d, time)
#'
#' Utime <- sort(unique(d$time))
#' 
#' ## fit cox models
#' m.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = d)
#' mStrata.CSC <- CSC(Hist(time, event) ~ strata(X1) + X2 + X3, data = d)
#'
#' ## compare models
#' outIF <- influenceTest(list(m.CSC, mStrata.CSC), 
#'              cause = 1, newdata = unique(d[,.(X1,X2,X3)]), times = Utime[1:5])
#' confint(outIF)

## * influenceTest.list (code)
#' @rdname influenceTest
#' @export
influenceTest.list <- function(object, newdata, times, type, cause,
                               keep.newdata = TRUE, keep.strata = FALSE, ...){
  
    validCox <- c("coxph","cph","phreg")
    name.model <- names(object)
    n.model <- length(object)
    if(is.null(name.model)){
       name.model <- paste0("model ",1:n.model)
    }
    ls.pred <- lapply(1:n.model, function(iM){
        iClass <- class(object[[iM]])
        
      if(any(validCox %in% iClass)){ # Cox
          pred <- predictCox(object[[iM]],
                             type = type,
                             newdata = newdata,
                             times = times,
                             iid = TRUE,
                             band = FALSE,
                             se = TRUE,
                             keep.newdata = keep.newdata,
                             keep.strata = keep.strata)
      }else if("CauseSpecificCox" %in% iClass){ # CSC
          pred <- predict(object[[iM]],
                          cause = cause,
                          newdata = newdata,
                          times = times,
                          iid = TRUE,
                          band = FALSE,
                          se = TRUE,
                          keep.newdata = keep.newdata,
                          keep.strata = keep.strata)
          pred$type <- "absRisk"
      }else{
          stop("can only handle Cox and Cause specific Cox models \n")
      }

      ## add call
      pred$name <- name.model[iM]
      pred$call <- object[[iM]]$call
      
      return(pred)
  })
    out <- influenceTest(ls.pred[[1]], ls.pred[[2]], ...)
    return(out)
}
                                        # }}}

## * influenceTest.default (code)
#' @rdname influenceTest
#' @export
influenceTest.default <- function(object,
                                  object2,
                                  band = TRUE,
                                  ...){
    
    ## ** check arguments
    if(class(object)=="predictCox"){
        if(object$type!=object2$type){
            stop("Cannot compare different types of prediction \n")
        }
        type <- object$type
        if(length(type)!=1){
            stop("Cannot analyse simulatenously several types of predictions \n",
                 "here: \"",paste0(type,collapse ="\" \""),"\"\n")
        }
    }else if(class(object)=="predictCSC"){
        type <- "absRisk"
    }else{
        stop("Can only deal with object of class \'predictCox\' or \'predictCSC\'")
    }
    estimator1 <- object[[type]]
    estimator2 <- object2[[type]]
    
    ## ** influence function
    type.iid <- paste0(type,".iid")
    if(type.iid %in% names(object) == FALSE){
        stop("No iid decomposition for type ",type," in object \n")
    }else{
        estimator1.iid <- object[[type.iid]]
    }
    if(type.iid %in% names(object2) == FALSE){
        stop("No iid decomposition for type ",type," in object2 \n")
    }else{
        estimator2.iid <- object2[[type.iid]]
    }
    
    ## check sample size
    n.obs <- length(estimator1.iid[,1,1])
    if(n.obs!=length(estimator2.iid[,1,1])){
        stop("Cannot compare estimates on different sample size \n")
    }
    
    ## time
    time <- object$time
    n.time <- length(time)

    ## compute the difference
    out <- list(name = c(object$name[1],object2$name[1]),
                call = list(object$call, object2$call),
                newdata = object$newdata,
                strata = object$strata,
                type = type,
                time = time,
                band = band)
    out$delta <- estimator1-estimator2
    out$delta.iid <- estimator1.iid-estimator2.iid
    out$delta.se <- t(sqrt(apply(out$delta.iid^2, MARGIN = 2:3, sum)))
    class(out) <- "influenceTest"

    if(band==FALSE){
        out$delta.iid <- NULL
    }
    
    return(out)
}

#----------------------------------------------------------------------
### influenceTest.R ends here
