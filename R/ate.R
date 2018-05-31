### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: maj 31 2018 (09:50) 
##           By: Brice Ozenne
##     Update #: 745
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Compute the average treatment effects via the g-formula
#'
#' @description Use the g-formula to estimate the average treatment
#'     effect based on Cox regression with or without competing risks
#' @param object outcome model which describes how event risk depends
#'     on treatment and covariates.  The object carry its own call and
#'     have a \code{predictRisk} method. See examples.
#' @param data data set in which to evaluate risk predictions based on
#'     the outcome model
#' @param formula For analyses with time-dependent covariates, the response formula. See examples.
#' @param treatment name of the treatment variable
#' @param contrasts the levels of the treatment variable to be
#'     compared
#' @param strata Strata variable on which to compute the average risk. Incompatible with treatment. Experimental.
#' @param times time points at which to evaluate risks
#' @param cause the cause of interest
#' @param landmark for models with time-dependent covariates the landmark time(s) of evaluation.
#'        In this case, argument \code{time} may only be one value and for the prediction of risks
#'        it is assumed that that the covariates do not change between landmark and landmark+time.
#' @param se Logical. If \code{TRUE} compute standard errors and confidence intervals
#' @param band Logical. If \code{TRUE} compute confidence bands across time points.
#' @param B the number of bootstrap replications used to compute the
#'     confidence intervals. If it equals 0, then Wald-type confidence
#'     intervals are computed.  They rely on the standard error
#'     estimated using the influence function of the estimator.
#' @param seed An integer used to generate seeds for bootstrap and to
#'     achieve reproducible results.
#' @param handler parallel handler for bootstrap. Either "mclapply" or
#'     "foreach". If "foreach" use \code{doParallel} to create a cluster.
#' @param mc.cores Passed to \code{parallel::mclapply} or
#'     \code{doParallel::registerDoParallel}. The number of cores to use, i.e. at
#'     most how many child processes will be run simultaneously.  The
#'     option is initialized from environment variable MC_CORES if
#'     set.
#' @param verbose Logical. If \code{TRUE} inform about estimated run
#'     time.
#' @param store.iid Implementation used to estimate the standard error. Can be \code{"full"} or \code{"minimal"}.
#' \code{"minimal"} requires less memory but can only estimate the standard for the difference between treatment effects (and not for the ratio).
#' @param ... passed to predictRisk
#'
#' @author Brice Ozenne \email{broz@@sund.ku.dk} and Thomas Alexander
#'     Gerds \email{tag@@biostat.ku.dk}
#' 
#' @examples 
#' library(survival)
#' library(rms)
#' 
#' set.seed(10)
#' 
#' #### Survival settings  ####
#' #### ATE with Cox model ####
#'
#' ## generate data
#' n <- 100
#' dtS <- sampleData(n, outcome="survival")
#' dtS$time <- round(dtS$time,1)
#' dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
#'
#' ## estimate the Cox model
#' fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' ## compute the ATE at times 5, 6, 7, and 8 using X1 as the treatment variable
#' \dontrun{
#' ## only punctual estimate (argument se = FALSE)
#' ateFit1a <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = FALSE)
#'
#' ## standard error / confidence intervals computed using the influence function
#' ## (argument se = TRUE and B = 0)
#' ateFit1b <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = TRUE, B = 0)
#'
#' ## bootstrap confidence intervals: studentized Wald type 
#' ateFit1c <- ate(fit, data = dtS, treatment = "X1", times = 5,
#'                seed=3,se = TRUE, B = 100)
#' ## bootstrap confidence intervals: studentized Wald type 
#' ateFit1d <- ate(fit, data = dtS, treatment = "X1", times = 5,
#'                 seed=3,bootci.method="quantile",se = TRUE, B = 100)
#'
#' ## same as before with in addition the confidence bands for the ATE
#' ## (argument band = TRUE)
#' ateFit1c <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = TRUE, band = TRUE, B = 0)
#' 
#' ## standard error / confidence intervals computed using 100 boostrap samples
#' ## (argument se = TRUE and B = 100) 
#' ateFit1d <- ate(fit, data = dtS, treatment = "X1",
#'                 times = 5:8, se = TRUE, B = 100)
#' ## NOTE: for real applications 100 bootstrap samples is not enougth 
#'
#' ## same but using 2 cpus for generating and analyzing the boostrap samples
#' ## (parallel computation, argument mc.cores = 2) 
#' ateFit1e <- ate(fit, data = dtS, treatment = "X1",
#'                 times = 5:8, se = TRUE, B = 100, mc.cores = 2)
#' }
#'
#' #### Survival settings without censoring ####
#' #### ATE with glm                        ####
#' 
#' ## generate data
#' n <- 100
#' dtS <- sampleData(n, outcome="survival")
#' dtS[, event5 := eventtime<=5]
#' dtS[, X2 := as.numeric(X2)]
#' 
#' ## estimate the Cox model
#' fit <- glm(formula = event5 ~ X1+X2, data=dtS, family = "binomial")
#' 
#' ## compute the ATE at times 5 using X1 as the treatment variable
#' ## only punctual estimate (argument se = FALSE)
#' ateFit1a <- ate(fit, data = dtS, treatment = "X1", times = 5,
#'                se = FALSE)
#' ateFit1a
#'
#' \dontrun{
#' ## standard error / confidence intervals computed using the influence function
#' ateFit1b <- ate(fit, data = dtS, treatment = "X1", times = 5,
#'                se = TRUE, B = 0)
#' ateFit1b
#'
#' ## standard error / confidence intervals computed using 100 boostrap samples
#' ateFit1d <- ate(fit, data = dtS, treatment = "X1",
#'                 times = 5, se = TRUE, B = 100)
#' ateFit1d
#' 
#' ## using lava
#' ateLava <- estimate(fit, function(p, data){
#' a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X2"] ;
#' R.X11 <- expit(a + b + c * data[["X2"]])
#' R.X10 <- expit(a + c * data[["X2"]])
#' list(risk0=R.X10,risk1=R.X11,riskdiff=R.X10-R.X11)},
#' average=TRUE)
#' ateLava
#' }
#' 
#' #### Competing risks settings               ####
#' #### ATE with cause specific Cox regression ####
#'
#' \dontrun{
#' ## generate data
#' n <- 500
#' dt <- sampleData(n, outcome="competing.risks")
#' dt$time <- round(dt$time,1)
#' dt$X1 <- factor(rbinom(n, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))
#'
#' ## estimate cause specific Cox model
#' fitCR <-  CSC(Hist(time,event)~ X1+X8,data=dt,cause=1)
#' 
#' ## compute the ATE at times 5, 6, 7, and 8 using X1 as the treatment variable
#' ateFit2a <- ate(fitCR, data = dt, treatment = "X1", times = 5:8, cause = 1,
#'                se = FALSE)
#'
#' ## standard error / confidence intervals computed using the influence function
#' ## (argument se = TRUE and B = 0)
#' ateFit2b <- ate(fitCR, data = dt, treatment = "X1", times = 5:8, cause = 1,
#'                se = TRUE, B = 0)
#'
#' ## same as before with in addition the confidence bands for the ATE
#' ## (argument band = TRUE)
#' ateFit2c <- ate(fitCR, data = dt, treatment = "X1", times = 5:8, cause = 1,
#'                se = TRUE, band = TRUE, B = 0)
#' 
#' ## standard error / confidence intervals computed using 100 boostrap samples
#' ## (argument se = TRUE and B = 100) 
#' ateFit2d <- ate(fitCR, data = dt, treatment = "X1", times = 5:8, cause = 1,
#'                 se = TRUE, B = 100)
#' ## NOTE: for real applications 100 bootstrap samples is not enougth 
#'
#' ## same but using 2 cpus for generating and analyzing the boostrap samples
#' ## (parallel computation, argument mc.cores = 2) 
#' ateFit2e <- ate(fitCR, data = dt, treatment = "X1", times = 5:8, cause = 1,
#'                 se = TRUE, B = 100, mc.cores = 2)
#' }
#' 
#' #### time-dependent covariates ###
#' \dontrun{
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ celltype+karno + age + trt, veteran)
#' vet2 <- survSplit(Surv(time, status) ~., veteran,
#'                        cut=c(60, 120), episode ="timegroup")
#' fitTD <- coxph(Surv(tstart, time, status) ~ celltype+karno + age + trt,
#'                data= vet2,x=1)
#' set.seed(16)
#' resVet <- ate(fitTD,formula=Hist(entry=tstart,time=time,event=status)~1,
#'           data = vet2, treatment = "celltype", contrasts = NULL,
#'         times=5,verbose=1,
#'         landmark = c(0,30,60,90), cause = 1, B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resVet
#' }
#' 
#' \dontrun{
#' set.seed(137)
#' d=sampleDataTD(127)
#' library(survival)
#' d[,status:=1*(event==1)]
#' ## ignore competing risks
#' cox1TD <- coxph(Surv(start,time, status,type="counting") ~ X3+X5+X6+X8, data=d)
#' resTD1 <- ate(cox1TD,formula=Hist(entry=start,time=time,event=status)~1,
#'         data = d, treatment = "X3", contrasts = NULL,
#'         times=.5,verbose=1,
#'         landmark = c(0,0.5,1), B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resTD1
#' ## adjust for competing risks
#' cscTD <- CSC(Hist(time=time, event=event,entry=start) ~ X3+X5+X6+X8, data=d)
#' set.seed(16)
#' resTD <- ate(cscTD,formula=Hist(entry=start,time=time,event=event)~1,
#'         data = d, treatment = "X3", contrasts = NULL,
#'         times=.5,verbose=1,
#'         landmark = c(0,0.5,1), cause = 1, B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resTD
#' }
#' @export
ate <- function(object,
                data,
                formula,
                treatment,
                strata = NULL,
                contrasts = NULL,
                times,
                cause,
                landmark,
                se = TRUE,
                iid = FALSE,
                band = FALSE,
                B = 0,
                seed,
                handler = "foreach",
                mc.cores = 1,
                verbose = TRUE,
                store.iid = "full",
                ...){
  
  meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b <- NULL
  .=.I <- NULL
  diff.se=ratio.se=.GRP=lower=upper=diff.lower=diff.upper=diff.p.value=ratio.lower=ratio.upper=ratio.p.value <- NULL
   
    handler <- match.arg(handler, c("foreach","mclapply","snow","parallel"))
                                        # {{{ checking for time-dependent covariates (left-truncation)
    TD <- switch(class(object)[[1]],"coxph"=(attr(object$y,"type")=="counting"),
                 "CauseSpecificCox"=(attr(object$models[[1]]$y,"type")=="counting"),FALSE)
    if (TD){
        if (missing(formula))
            stop("Need formula to do landmark analysis.")
        if (missing(landmark))
            stop("Need landmark time(s) to do landmark analysis.")
        if(length(times)!=1){
            stop("In settings with time-dependent covariates argument 'time' must be a single value, argument 'landmark' may be a vector of time points.")
        }
        Gformula <- Gformula_TD
    }else{
        landmark <- NULL
        Gformula <- Gformula_TI
    }
                                        # }}}
                                        # {{{ Prepare
    dots <- list(...)
    if(se==0 && B>0){
        warning("argument 'se=0' means 'no standard errors' so number of bootstrap repetitions is forced to B=0.")
    }
    if(band && B>0){
        stop("Confidence bands cannot be computed when using the bootstrap approach \n",
             "Either set argument \'band\' to FALSE to not compute the confidence bands \n",
             "or set argument \'B\' to 0 to use the estimate of the asymptotic distribution instead of the bootstrap\n")
    }
    if(!is.null(treatment)){
        if(length(treatment) != 1){
            stop("Argument treatment should have length 1. \n")
        }
        if(treatment %in% names(data) == FALSE){
            stop("The data set does not seem to have a variable ",treatment," (argument: treatment). \n")
        }
        if(!is.null(strata) ){
            stop("Argument strata must be NULL when argument treatment is specified. \n")
        }
        strata <- treatment
    }else{
        if(is.null(strata) ){
            stop("Argument strata must refer to a variable in the data when argument treatment is NULL. \n")

        }
        if(length(strata) != 1){
            stop("Argument strata should have length 1. \n")
        }
        if(strata %in% names(data) == FALSE){
            stop("The data set does not seem to have a variable ",strata," (argument: strata). \n")
        }
        if(B > 0){
            stop("Boostrap resampling is not available when argument strata is specified. \n")
        }
        if(TD){
            stop("Landmark analysis is not available when argument strata is specified. \n")
        }
    }
    test.CR <- !missing(cause) # test whether the argument cause has been specified, i.e. it is a competing risk model
    if(test.CR==FALSE){cause <- NA}
  
    if(B==0 && (se || band)){
        validClass <- c("CauseSpecificCox","coxph","cph","phreg","glm")
        if(all(validClass %in% class(object) == FALSE)){
            stop("Standard error based on the influence function only implemented for \n",
                 paste(validClass, collapse = " ")," objects \n",
                 "set argument \'B\' to a positive integer to use a boostrap instead \n")
        }
    }
    if(!is.null(treatment)){
        data[[treatment]] <- factor(data[[treatment]])
    
        if(is.null(contrasts)){
            levels <- levels(data[[treatment]])
            contrasts <- levels(data[[treatment]])
            ## if (length(contrasts)>50) warning("Treatment variable has more than 50 levels.\nIf this is not a mistake,
            ## you should use the argument `contrasts'.")
        }else{levels <- contrasts}
    }else{
        data[[strata]] <- factor(data[[strata]])        
        levels <- levels(data[[strata]])
        contrasts <- levels(data[[strata]])
    }
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    n.obs <- NROW(data)
  
  # }}}
  
    # {{{ Checking the model
  
    ## for predictRisk S3-method
    allmethods <- utils::methods(predictRisk)
    candidateMethods <- paste("predictRisk",class(object),sep=".")
    if (all(match(candidateMethods,allmethods,nomatch=0)==0))
        stop(paste("Could not find predictRisk S3-method for ",class(object),collapse=" ,"),sep="")
    ## for compatibility with resampling
    if(is.null(object$call))
        stop(paste("The object does not contain its own call, which is needed to refit the model in the bootstrap loop."))
    # }}}

    # {{{ Point estimate
    estimateTime <- system.time(
        pointEstimate <- Gformula(object=object,
                                  data=data,
                                  treatment=treatment,
                                  strata=strata,
                                  contrasts=contrasts,
                                  times=times,
                                  cause=cause,
                                  landmark=landmark,
                                  n.contrasts = n.contrasts,
                                  levels = levels,
                                  dots)
    )
    # }}}

    # {{{ Confidence interval    

    if(se || band){
        if (TD){
            key1 <- c("Treatment","landmark")
            key2 <- c("Treatment.A","Treatment.B","landmark")
        }
        else{
            key1 <- c("Treatment","time")
            key2 <- c("Treatment.A","Treatment.B","time")
        }
    
        if(B>0){
                                        # {{{ Bootstrap
            if (verbose==TRUE){ ## display
                message(paste0("Approximated bootstrap netto run time (without time for copying data to cores):\n",
                               round(estimateTime["user.self"],2),
                               " seconds times ",
                               B,
                               " bootstraps / ",
                               mc.cores,
                               " cores = ",
                               round(estimateTime["user.self"]*B/mc.cores,2)," seconds.\n"))
            }

            vec.pointEstimate <- c(pointEstimate$meanRisk$meanRisk,
                                   pointEstimate$riskComparison$diff,
                                   pointEstimate$riskComparison$ratio)
            names(vec.pointEstimate) <- c(pointEstimate$meanRisk[,paste0("meanRisk:",Treatment,":",time)],
                                          pointEstimate$riskComparison[,paste0("compRisk:diff:",Treatment.A,":",Treatment.B,":",time)],
                                          pointEstimate$riskComparison[,paste0("compRisk:ratio:",Treatment.A,":",Treatment.B,":",time)]
                                          )

            resBoot <- calcBootATE(object,
                                   pointEstimate = vec.pointEstimate,
                                   Gformula = Gformula,
                                   data = data,
                                   treatment = treatment,
                                   contrasts = contrasts,
                                   times = times,
                                   cause = cause,
                                   landmark = landmark,
                                   n.contrasts = n.contrasts,
                                   levels = levels,
                                   dots = dots,
                                   n.obs = n.obs,
                                   handler = handler,
                                   B = B,
                                   seed = seed,
                                   mc.cores = mc.cores,
                                   verbose = verbose)

            bootseeds <- resBoot$bootseeds
            resBoot <- resBoot$boot
            outSE <- NULL   
                                        # }}}
        } else {
                                        # {{{ compute standard error and quantiles via the influence function

            if(!is.null(landmark)){
                stop("Calculation of the standard errors via the influence function not implemented for time dependent covariates \n")
            }
            outSE <- calcSeATE(object,
                               data = data,
                               times = times,
                               cause = cause,
                               treatment = treatment,
                               contrasts = contrasts,
                               strata = strata,
                               n.contrasts = n.contrasts,
                               levels = levels,
                               n.times = n.times,
                               n.obs = n.obs,
                               pointEstimate = pointEstimate,
                               se = se,
                               iid = (band+iid)>0, 
                               store.iid = store.iid)

            bootseeds <- NULL
            resBoot <- NULL          
   
                                        # }}}
        }
    } else{
        bootseeds <- NULL
        resBoot <- NULL
        outSE <- NULL

                                        # }}}
    }
                                        # {{{ output object
    if(is.null(treatment)){
        setnames(pointEstimate$meanRisk,
                 old = "Treatment",
                 new = strata)
        setnames(pointEstimate$riskComparison,
                 old = c("Treatment.A","Treatment.B"),
                 new = paste0(strata,c(".A",".B")))
    }
    
    out <- c(list(meanRisk = pointEstimate$meanRisk,
                  riskComparison = pointEstimate$riskComparison),
             outSE,
             list(treatment = treatment,
                  contrasts = contrasts,
                  times = times,
                  se = se,
                  TD = TD,
                  B = B,
                  band = band,
                  boot = resBoot,
                  seeds = bootseeds)
             )
  
    class(out) <- c("ate",class(object))
    return(out)
                                        # }}}

}

# {{{ Gformula: time dependent covariates
Gformula_TD <- function(object,
                        data,
                        treatment,
                        contrasts,
                        times,
                        landmark,
                        cause,
                        n.contrasts,
                        levels,
                        ...){

    Treatment <- Treatment.B <- meanRisk <- ratio <- NULL ## [:forCRANcheck:]
    
    response <- eval(formula[[2]],envir=data)
    time <- response[,"time"]
    entry <- response[,"entry"]
    if(class(object)[[1]]=="coxph"){
        riskhandler <- "predictRisk.coxphTD"
    }else{
        riskhandler <- "predictRisk.CSCTD"
    }
    ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
    dt.meanRisk <- data.table::rbindlist(lapply(1:n.contrasts,function(i){
        data.i <- data
        data.i[[treatment]] <- factor(contrasts[i], levels = levels)
        data.table::rbindlist(lapply(landmark,function(lm){
            atrisk <- (entry <= lm & time >= lm)
            risk.i <- colMeans(do.call(riskhandler,
                                       args = list(object,
                                                   newdata = data.i[atrisk,],
                                                   times = times,
                                                   cause = cause,
                                                   landmark=lm,
                                                   ...)))
            data.table::data.table(Treatment=contrasts[[i]],time=times,landmark=lm,meanRisk=risk.i)
        }))
    }))
    riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){
            ## compute differences between all pairs of treatments
            RC <- dt.meanRisk[Treatment==contrasts[[i]]]
            setnames(RC,"Treatment","Treatment.A")
            RC[,Treatment.B:=contrasts[[j]]]
            RC[,diff:=meanRisk-dt.meanRisk[Treatment==contrasts[[j]],meanRisk]]
            RC[,ratio:=meanRisk/dt.meanRisk[Treatment==contrasts[[j]],meanRisk]]
            RC[,meanRisk:=NULL]
            RC[]
        }))}))
    out <- list(meanRisk = dt.meanRisk,
                riskComparison = riskComparison,
                treatment = treatment,
                strata = strata)
    return(out)
}

# }}}

# {{{ Gformula: time independent covariates
Gformula_TI <- function(object,
                        data,
                        treatment,
                        strata,
                        contrasts,
                        times,
                        landmark,
                        cause,
                        n.contrasts,
                        levels,
                        ...){

    meanRisk <- lapply(1:n.contrasts,function(i){ ## i <- 1
        ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
        if(!is.null(treatment)){
            data.i <- data
            data.i[[treatment]] <- factor(contrasts[i], levels = levels)
        }else{
            data.i <- data[data[[strata]]==contrasts[i]]
        }
        allrisks <- do.call("predictRisk",
                            args = list(object, newdata = data.i, times = times, cause = cause,...))
        if(!is.matrix(allrisks)){allrisks <- cbind(allrisks)} 
        risk.i <- colMeans(allrisks)
        risk.i
    })

        riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){ ## i <- 1
            data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){ ## j <- 2
                ## compute differences between all pairs of treatments
                data.table(Treatment.A=contrasts[[i]],
                           Treatment.B=contrasts[[j]],
                           time = times,
                           diff=meanRisk[[i]]-meanRisk[[j]],
                           ratio=meanRisk[[i]]/meanRisk[[j]])
            }))}))

    ## reshape for export
    name.strata <- unlist(lapply(1:n.contrasts, function(c){rep(contrasts[c],length(meanRisk[[c]]))}))

    out <- list(meanRisk = data.table(Treatment=name.strata,
                                      time = times,
                                      meanRisk=unlist(meanRisk)),
                riskComparison = riskComparison)
    return(out)            
}
# }}}

            

#----------------------------------------------------------------------
### ate.R ends here
