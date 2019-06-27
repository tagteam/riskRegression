### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: jun 27 2019 (13:49) 
##           By: Brice Ozenne
##     Update #: 945
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * ate (documentation)
#' @title Compute the Average Treatment Effects Via the G-formula
#' @description Use the g-formula to estimate the average treatment
#'     effect based on Cox regression with or without competing risks.
#' @name ate
#' 
#' @param object.event Outcome model which describes how event risk depends
#'     on treatment and covariates. The object carry its own call and
#'     have a \code{predictRisk} method. See examples.
#' @param object.treatment Treatment model which describes how treatment depends
#'     on covariates. The object must be a \code{glm} object (logistic regression). See examples.
#' @param object.censor Censoring model which describes how censoring depends
#'     on treatment and covariates. The object must be a \code{coxph} or \code{cph} object. See examples.
#' @param data [data.frame or data.table] Data set in which to evaluate risk predictions
#' based on the outcome model
#' 
#' @param formula For analyses with time-dependent covariates, the response formula. See examples.
#' @param treatment [character] Name of the treatment variable.
#' @param contrasts [character] The levels of the treatment variable to be compared.
#' @param strata [character] Strata variable on which to compute the average risk.
#' Incompatible with treatment. Experimental.
#' 
#' @param times [numeric vector] Time points at which to evaluate average treatment effects.
#' @param cause [integer/character] the cause of interest.
#' @param landmark for models with time-dependent covariates the landmark time(s) of evaluation.
#'        In this case, argument \code{time} may only be one value and for the prediction of risks
#'        it is assumed that that the covariates do not change between landmark and landmark+time.
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param B [integer, >0] the number of bootstrap replications used to compute the confidence intervals.
#' If it equals 0, then the influence function is used to compute Wald-type confidence intervals/bands.
#' @param seed [integer, >0] sed number used to generate seeds for bootstrap
#' and to achieve reproducible results.
#' @param handler [character] Parallel handler for bootstrap.
#' either \code{"foreach"}, \code{"mclapply"}, \code{"snow"} or \code{"multicore"}.
#' if "foreach" use \code{doparallel} to create a cluster.
#' @param mc.cores [integer, >0] The number of cores to use,
#' i.e., the upper limit for the number of child processes that run simultaneously.
#' Passed to \code{parallel::mclapply} or \code{doparallel::registerdoparallel}.
#' The option is initialized from environment variable mc_cores if set.
#' @param cl A parallel socket cluster used to perform cluster calculation in parallel.
#' Output by \code{parallel::makeCluster}.
#' The packages necessary to run the computations (e.g. riskRegression) must already be loaded on each worker.
#' @param verbose [logical] If \code{TRUE} inform about estimated run time.
#' @param store.iid [character] Implementation used to estimate the standard error.
#' Can be \code{"full"} or \code{"minimal"}.
#' \code{"minimal"} requires less memory but can only estimate the standard for the difference between treatment effects (and not for the ratio).
#' @param ... passed to predictRisk
#'
#' @author Brice Ozenne \email{broz@@sund.ku.dk}
#' and Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#'
#' @seealso
#' \code{\link{confint.ate}} to compute confidence intervals/bands.
#' \code{\link{autoplot.ate}} to display the average risk.
#' \code{\link{ateRobust}} to make the estimator doubly robust 
## * ate (examples)
#' @rdname ate
#' @examples 
#' library(survival)
#' library(rms)
#' library(prodlim)
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
#' ## only point estimate (argument se = FALSE)
#' ateFit1a <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = TRUE)
#'
#' ## standard error / confidence intervals computed using the influence function
#' ## (argument se = TRUE and B = 0)
#' ateFit1b <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = TRUE, B = 0)
#' 
#' ## same as before with in addition the confidence bands for the ATE
#' ## (argument band = TRUE)
#' ateFit1c <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                se = TRUE, band = TRUE, B = 0)
#'
#' ## bootstrap confidence intervals
#' ateFit1c <- ate(fit, data = dtS, treatment = "X1", times = 5,
#'                seed = 3, se = TRUE, B = 100)
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
#' dtB <- sampleData(n, outcome="binary")
#' dtB[, X2 := as.numeric(X2)]
#' 
#' ## estimate a logistic regression model
#' fit <- glm(formula = Y ~ X1+X2, data=dtB, family = "binomial")
#' 
#' ## compute the ATE using X1 as the treatment variable
#' ## only point estimate (argument se = FALSE)
#' ateFit1a <- ate(fit, data = dtB, treatment = "X1", se = FALSE)
#' ateFit1a
#'
#' \dontrun{
#' ## standard error / confidence intervals computed using the influence function
#' ateFit1b <- ate(fit, data = dtB, treatment = "X1",
#'                times = 5, ## just for having a nice output not used in computations
#'                se = TRUE, B = 0)
#' ateFit1b
#'
#' ## standard error / confidence intervals computed using 100 boostrap samples
#' ateFit1d <- ate(fit, data = dtB, treatment = "X1",
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
#'         landmark = c(0,30,60,90), cause = 1, B = 4, se = 1,
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

## * ate (code)
#' @rdname ate
#' @export
ate <- function(object.event,
                object.censor,
                object.treatment,
                data,
                formula,
                treatment,
                strata = NULL,
                contrasts = NULL,
                times,
                cause = NA,
                landmark,
                se = TRUE,
                iid = FALSE,
                band = FALSE,
                B = 0,
                confint = (se+band)>0,
                seed,
                handler = "foreach",
                mc.cores = 1,
                cl = NULL,
                verbose = TRUE,
                store.iid = "full",
                ...){

    dots <- list(...)

    ## for CRAN tests
    meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b <- NULL
    .=.I <- NULL
    diff.se=ratio.se=.GRP=lower=upper=diff.lower=diff.upper=diff.p.value=ratio.lower=ratio.upper=ratio.p.value <- NULL

    ## ** initialize arguments
    init <- ate_initArgs(object.event = object.event,
                         treatment = treatment,
                         strata = strata,
                         times = times,
                         handler = handler)
    times <- init$times
    handler <- init$handler
    strata <- init$strata
    formula <- init$formula
    landmark <- init$landmark
    TD <- init$TD
    fct.pointEstimate <- init$fct.pointEstimate
    n.train <- init$n.train
    
    ## ** check consistency of the user arguments
    check <- ate_checkArgs(object.event = object.event,
                           object.censor = object.censor,
                           object.treatment = object.treatment,
                           data = data,
                           formula = formula,
                           treatment = treatment,
                           strata = strata,
                           contrasts = contrasts,
                           times = times,
                           cause = cause,
                           landmark = landmark,
                           se = se,
                           iid = iid,
                           band = band,
                           B = B,
                           confint = confint,
                           seed = seed,
                           handler = handler,
                           mc.cores = mc.cores,
                           cl = cl,
                           verbose = verbose,
                           store.iid = store.iid,
                           augment.cens = augment.cens,
                           TD = TD,
                           n.train = n.train)

    ## ** Prepare
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
  
    ## ** Point estimate
    args.pointEstimate <- list(object.event=object.event,
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
    if (TD){
        args.pointEstimate <- c(args.pointEstimate,list(formula=formula))
    }
    estimateTime <- system.time(pointEstimate <- do.call(fct.pointEstimate, args.pointEstimate))

    ## ** Confidence intervals
    if(se || band || iid){
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
            ### *** Bootstrap
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

            resBoot <- calcBootATE(object.event,
                                   pointEstimate = vec.pointEstimate,
                                   fct.pointEstimate = fct.pointEstimate,
                                   data = data,
                                   formula=formula,
                                   TD=TD,
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
                                   cl = cl,
                                   verbose = verbose)

            bootseeds <- resBoot$bootseeds
            resBoot <- resBoot$boot
            outSE <- NULL   
                                        # }}}
        } else {
                                        # {{{ compute standard error and quantiles via the influence function
            ## *** Delta method
            outSE <- calcSeATE(object.event,
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
                               export = c("iid"[(iid+band)>0],"se"[(se+band)>0]),
                               store.iid = store.iid)

            pointEstimate$meanRisk[["meanRisk.se"]] <- outSE$meanRisk.se[,1]
            pointEstimate$riskComparison[["diff.se"]] <- outSE$diffRisk.se[,1]
            pointEstimate$riskComparison[["ratio.se"]] <- outSE$ratioRisk.se[,1]            
            data.table::setcolorder(pointEstimate$riskComparison, neworder = c(names(pointEstimate$riskComparison)[1:3],
                                                                               "diff","diff.se","ratio","ratio.se"))
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
    ## ** export
    if(is.null(treatment)){
        setnames(pointEstimate$meanRisk,
                 old = "Treatment",
                 new = strata)
        setnames(pointEstimate$riskComparison,
                 old = c("Treatment.A","Treatment.B"),
                 new = paste0(strata,c(".A",".B")))
    }
    
    out <- list(meanRisk = pointEstimate$meanRisk,
                riskComparison = pointEstimate$riskComparison,
                meanRisk.iid = outSE$meanRisk.iid,
                diffRisk.iid = outSE$diffRisk.iid,
                ratioRisk.iid = outSE$ratioRisk.iid,
                treatment = treatment,
                contrasts = contrasts,
                times = times,
                se = se,
                TD = TD,
                B = B,
                band = band,
                boot = resBoot,
                seeds = bootseeds)

  
    class(out) <- c("ate",class(object.event))
    if(confint){
        out <- stats::confint(out)
    }
    if(band[[1]] && se[[1]]==FALSE){
        out$meanRisk[["meanRisk.se"]] <- NULL
        out$riskComparison[["diff.se"]] <- NULL
        out$riskComparison[["ratio.se"]] <- NULL
    }
    if(band[[1]] && iid[[1]]==FALSE){
        out[paste0(c("mean","diff","ratio"),"Risk.iid")] <- NULL
    }

    return(out)

}

## * ate_initArgs
ate_initArgs <- function(object.event,
                         treatment,
                         strata,
                         times,
                         handler){

    ## ** user-defined arguments
    ## times
    if(inherits(object.event,"glm") && missing(times)){
        times <- NA
    }

    ## handler
    handler <- match.arg(handler, c("foreach","mclapply","snow","multicore"))

    ## strata
    if(!is.null(treatment) & is.null(strata)){
        strata <- treatment
    }

    ## ** deduced from user defined arguments
    ## presence of time dependent covariates
    TD <- switch(class(object.event)[[1]],
                 "coxph"=(attr(object.event$y,"type")=="counting"),
                 "CauseSpecificCox"=(attr(object.event$models[[1]]$y,"type")=="counting"),
                 FALSE)
    if(TD){
        fct.pointEstimate <- Gformula_TD
    }else{
        landmark <- NULL
        formula <- NULL
        fct.pointEstimate <- Gformula_TI
    }

    n.train <- coxN(object.event)
    
    ## ** output
    return(list(times = times,
                handler = handler,
                strata = strata,
                formula = formula,
                landmark = landmark,
                TD = TD,
                fct.pointEstimate = fct.pointEstimate,
                n.train = n.train))
}

## * ate_checkArgs
ate_checkArgs <- function(object.event,
                          object.censor,
                          object.treatment,
                          data,
                          formula,
                          treatment,
                          strata,
                          contrasts,
                          times,
                          cause,
                          landmark,
                          se,
                          iid,
                          band,
                          B,
                          confint,
                          seed,
                          handler,
                          mc.cores,
                          cl,
                          verbose,
                          store.iid,
                          augment.cens,
                          TD,
                          n.train){

    ## ** times
    if(inherits(object.event,"glm") && length(times)!=1){
        warning("Argument \'times\' has no effect when using a glm object \n",
                "It should be set to NA \n")        
    }

    ## ** object.event
    ## is there a predict method?
    allmethods <- utils::methods(predictRisk)
    candidateMethods <- paste("predictRisk",class(object.event),sep=".")
    if (all(match(candidateMethods,allmethods,nomatch=0)==0)){
        stop(paste("Could not find predictRisk S3-method for ",class(object.event),collapse=" ,"),sep="")
    }
    
    ## has the delta method been implemented?
    if(B==0 && (se || band || iid)){
        validClass <- c("CauseSpecificCox","coxph","cph","phreg","glm")
        if(all(validClass %in% class(object.event) == FALSE)){
            stop("Standard error based on the influence function only implemented for \n",
                 paste(validClass, collapse = " ")," objects \n",
                 "set argument \'B\' to a positive integer to use a boostrap instead \n")
        }
        if(n.train[1]!=NROW(data)){
            stop("Argument \'data\' must contain the dataset used to fit the object (when \'se\', \'band\', or \'iid\' is TRUE)\n")
        }
    }
   
    ## ** bootstrap
    if(B>0){
        if(se==FALSE){
            warning("Argument 'se=0' means 'no standard errors' so number of bootstrap repetitions is forced to B=0.")
        }
        if(iid==TRUE){
            stop("Influence function cannot be computed when using the bootstrap approach \n",
                 "Either set argument \'iid\' to FALSE to not compute the influence function \n",
                 "or set argument \'B\' to 0 \n")
        }
        if(band==TRUE){
            stop("Confidence bands cannot be computed when using the bootstrap approach \n",
                 "Either set argument \'band\' to FALSE to not compute the confidence bands \n",
                 "or set argument \'B\' to 0 to use the estimate of the asymptotic distribution instead of the bootstrap\n")
        }
        if(is.null(object.event$call)){
            stop("The object does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        ## if((Sys.info()["sysname"] == "Windows") && (handler == "mclapply") && (mc.cores>1) ){
            ## stop("mclapply cannot perform parallel computations on Windows \n",
                 ## "consider setting argument handler to \"foreach\" \n")
        ## }
        max.cores <- parallel::detectCores()
        if(mc.cores > max.cores){
            stop("Not enough available cores \n","available: ",max.cores," | requested: ",mc.cores,"\n")
        }
    }

    ## ** delta method
    if(B==0 && (se|band|iid) && !is.null(landmark)){
        stop("Calculation of the standard errors via the influence function not implemented for time dependent covariates \n")
    }

    
    ## ** treatment & strata
    if(identical(treatment,strata)){
        if(length(treatment) != 1){
            stop("Argument treatment should have length 1. \n")
        }
        if(treatment %in% names(data) == FALSE){
            stop("The data set does not seem to have a variable ",treatment," (argument: treatment). \n")
        }
        if(is.numeric(data[[treatment]])){
            stop("The treatment variable must be a factor variable. \n",
                 "Convert treatment to factor, re-fit the object using this new variable and then call ate. \n")
        }
    }else{
        if(!is.null(treatment)){
            stop("Argument treatment must be NULL when strata is specified. \n") ## the case strata==treatment is hidden to the user
        }
        if(any(strata %in% names(data) == FALSE)){
            stop("The data set does not seem to have a variable \"",paste0(strata, collapse = "\" \""),"\" (argument: strata). \n")
        }
        if(length(strata) != 1){
            stop("Argument strata should have length 1. \n")
        }
        if(B > 0){
            stop("Boostrap resampling is not available when argument strata is specified. \n")
        }
        if(TD == TRUE){
            stop("Landmark analysis is not available when argument strata is specified. \n")
        }
    }
    
    ## ** time dependent covariances
    if (TD){
        if (missing(formula))
            stop("Need formula to do landmark analysis.")
        if (missing(landmark))
            stop("Need landmark time(s) to do landmark analysis.")
        if(length(times)!=1){
            stop("In settings with time-dependent covariates argument 'time' must be a single value, argument 'landmark' may be a vector of time points.")
        }
    }

    ## ** n.train
    ## when using CSC model, check that the same dataset has been used to train each Cox model
    ## in fact only check that via the size of the dataset
    if(length(unique(n.train))>1){
        stop("The same dataset must have been used to fit all models.")
    }
    
    ## ** output
    return(TRUE)
}


            

#----------------------------------------------------------------------
### ate.R ends here
