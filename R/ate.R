### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: sep 16 2019 (18:08) 
##           By: Brice Ozenne
##     Update #: 1347
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
#' @param event Outcome model which describes how event risk depends
#'     on treatment and covariates. The object carry its own call and
#'     have a \code{predictRisk} method. See examples.
#' @param treatment Treatment model which describes how treatment depends
#'     on covariates. The object must be a \code{glm} object (logistic regression) or the name of the treatment variable.
#'     See examples.
#' @param censor Censoring model which describes how censoring depends
#'     on treatment and covariates. The object must be a \code{coxph} or \code{cph} object. See examples.
#' @param data [data.frame or data.table] Data set in which to evaluate risk predictions
#' based on the outcome model
#' 
#' @param formula For analyses with time-dependent covariates, the response formula. See examples.
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
#'                se = FALSE)
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
#'
#' \dontrun{
#' ## standard error / confidence intervals computed using the influence function
#' ateFit1b <- ate(fit, data = dtB, treatment = "X1",
#'                times = 5, ## just for having a nice output not used in computations
#'                se = TRUE, B = 0)
#'
#' ## standard error / confidence intervals computed using 100 boostrap samples
#' ateFit1d <- ate(fit, data = dtB, treatment = "X1",
#'                 times = 5, se = TRUE, B = 100)
#' 
#' ## using the lava package
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
#' d[,X3:=as.factor(X3)]
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
ate <- function(event,
                treatment,
                censor = NULL,
                data,
                formula,
                strata = NULL,
                contrasts = NULL,
                times,
                cause = NA,
                landmark,
                se = TRUE,
                iid = FALSE,
                band = FALSE,
                B = 0,
                seed,
                handler = "foreach",
                mc.cores = 1,
                cl = NULL,
                verbose = 2,
                ...){

    dots <- list(...)

    ## for CRAN tests
    meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b <- NULL
    .=.I <- NULL
    diff.se=ratio.se=.GRP=lower=upper=diff.lower=diff.upper=diff.p.value=ratio.lower=ratio.upper=ratio.p.value <- NULL

    ## ** initialize arguments
    data.table::setDT(data)

    init <- ate_initArgs(object.event = event,
                         object.treatment = treatment,
                         object.censor = censor,
                         formula = formula,
                         landmark = landmark,
                         data = data,
                         times = times,
                         cause = cause,
                         handler = handler,
                         product.limit = dots$product.limit)

    object.event <- init$object.event
    object.treatment <- init$object.treatment
    object.censor <- init$object.censor
    times <- init$times
    handler <- init$handler
    formula <- init$formula
    landmark <- init$landmark
    treatment <- init$treatment
    cause <- init$cause
    TD <- init$TD
    estimator <- init$estimator
    fct.pointEstimate <- init$fct.pointEstimate
    n.obs <- init$n.obs
    n.censor <- init$n.censor
    eventVar.time <- init$eventVar.time
    eventVar.status <- init$eventVar.status
    level.censoring <- init$level.censoring
    censorVar.time <- init$censorVar.time
    censorVar.status <- init$censorVar.status
    dots$product.limit <- init$product.limit

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
                           seed = seed,
                           handler = handler,
                           mc.cores = mc.cores,
                           cl = cl,
                           verbose = verbose,
                           TD = TD,
                           n.censor = n.censor,
                           level.censoring = level.censoring,
                           estimator = estimator,
                           eventVar.time = eventVar.time,
                           eventVar.status = eventVar.status,
                           censorVar.time = censorVar.time,
                           censorVar.status = censorVar.status
                           )

    ## ** Prepare
    if(!is.null(strata)){
        data[[strata]] <- factor(data[[strata]])        
        if(is.null(contrasts)){
            contrasts <- levels(data[[strata]])
        }
        levels <- levels(data[[strata]]) ## necessary to get the correct factor format in case that not all levels are in contrast
    }else{
        data[[treatment]] <- factor(data[[treatment]])
    
        if(is.null(contrasts)){
            contrasts <- levels(data[[treatment]])
            ## if (length(contrasts)>50) warning("Treatment variable has more than 50 levels.\nIf this is not a mistake,
            ## you should use the argument `contrasts'.")
        }
        levels <- levels(data[[treatment]]) ## necessary to get the correct factor format in case that not all levels are in contrast
    }

    ## ** Overview of the calculations
    if(verbose>1){
        name.estimator <- switch(estimator,
                                 "Gformula" = "the G-formula",
                                 "GformulaTD" = "the G-formula with time dependent covariates",
                                 "IPTW" = "the IPTW estimator",
                                 "IPTW,IPCW" = "the IPTW, IPCW estimator",
                                 "AIPTW" = "the AIPTW estimator",
                                 "AIPTW,AIPCW" = "the AIPTW, AIPCW estimator"
                                 )

        if(!is.null(strata)){
            cat("Calculation of the average risk by strata \n", sep = "")
            cat(" - Estimator      : ",estimator,"\n",sep="")
            cat(" - Strata variable: ",strata,"\n",sep="")
        }else{
            cat("     Calculation of the average treatment effect \n\n", sep = "")
            cat(" - Estimator         : ",estimator,"\n",sep="")
            cat(" - Treatment variable: ",treatment," (",length(contrasts)," levels)\n",sep="")
            cat(" - Event variable    : ",eventVar.status," (cause: ",cause,", censoring: ",level.censoring,")\n",sep="")
            cat(" - Time variable     : ",eventVar.time,"\n",sep="")
        }
        cat("\n")
        cat(" - Point estimation:")
        
    }

    ## ** Point estimate
    args.pointEstimate <- c(list(object.event = object.event,
                               object.censor = object.censor,
                               object.treatment = object.treatment,
                               data=data,
                               treatment=treatment,
                               strata=strata,
                               contrasts=contrasts,
                               levels=levels,
                               times=times,
                               cause=cause,
                               landmark=landmark,
                               n.censor = n.censor,
                               level.censoring = level.censoring,
                               estimator = estimator,
                               eventVar.time = eventVar.time,
                               eventVar.status = eventVar.status,
                               censorVar.time = censorVar.time,
                               censorVar.status = censorVar.status,
                               return.iid = (se || band || iid) && (B==0)),
                               dots)
    if (TD){       
        args.pointEstimate <- c(args.pointEstimate,list(formula=formula))
    }
    ## note: system.time() seems to slow down the execution of the function, this is why Sys.time is used instead.
    if(B>0){
        tps1 <- Sys.time()
    }

    pointEstimate <- do.call(fct.pointEstimate, args.pointEstimate)
    if(B>0){
        tps2 <- Sys.time()
        estimateTime <- as.numeric(tps2-tps1)
    }
    
    if(verbose>1){cat(" done \n")}

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
            ## *** Bootstrap
            ## display start
            if (verbose>1){ 
                cat(" - Non-parametric bootstrap using ",B," samples and ",mc.cores," core",if(mc.cores>1){"s"},"\n", sep ="")
                cat("                            (expected time: ",round(estimateTime*B/mc.cores,2)," seconds)\n", sep = "")
            }

            ## prepare arguments
            vec.pointEstimate <- c(pointEstimate$meanRisk$meanRisk,
                                   pointEstimate$riskComparison$diff,
                                   pointEstimate$riskComparison$ratio)
            names(vec.pointEstimate) <- c(pointEstimate$meanRisk[,paste0("meanRisk:",.SD$Treatment,":",.SD$time)],
                                          pointEstimate$riskComparison[,paste0("compRisk:diff:",.SD$Treatment.A,":",.SD$Treatment.B,":",.SD$time)],
                                          pointEstimate$riskComparison[,paste0("compRisk:ratio:",.SD$Treatment.A,":",.SD$Treatment.B,":",.SD$time)]
                                          )

            ## run
            resBoot <- calcBootATE(args = args.pointEstimate,
                                   name.estimate = names(vec.pointEstimate),
                                   n.obs = n.obs,
                                   fct.pointEstimate = fct.pointEstimate,
                                   handler = handler,
                                   B = B,
                                   seed = seed,
                                   mc.cores = mc.cores,
                                   cl = cl,
                                   verbose = verbose)

            ## store
            boot.object <- list(t0 = vec.pointEstimate,
                                t = resBoot$boot,
                                R = B,
                                data = data,
                                seed = resBoot$bootseeds,
                                statistic = NULL,
                                sim = "ordinary",
                                call = quote(boot(data = XX, statistic = XX, R = XX)),
                                stype = "i",
                                strata = rep(1,n.obs),
                                weights = rep(1/n.obs,n.obs),
                                pred.i = NULL,  ## Omitted if m is 0 or sim is not "ordinary" (from doc of boot::boot)
                                L = NULL, ## only used when sim is "antithetic" (from doc of boot::boot)
                                ran.gen = NULL, ## only used when sim is "parametric" (from doc of boot::boot)
                                mle = NULL ## only used when sim is "parametric" (from doc of boot::boot)
                                )
            class(boot.object) <- "boot"
            bootseeds <- resBoot$bootseeds
            outIID <- NULL

            ## display end
            if (verbose>1){ 
                cat("\n")
            }

                                        # }}}
        } else {
                                        # {{{ compute standard error and quantiles via the influence function
            ## *** Delta method
            ## display start
            if (verbose>1){ 
                cat(" - Functional delta method: ")
            }

            if(is.null(attr(iid,"nuisance")) || (attr(iid,"nuisance")==TRUE)){ ## compute iid decomposition relative to the nuisance parameters
                

                ## add pre-computed quantities
                name.attributes <- setdiff(names(attributes(pointEstimate)),"names")
                for(iAttr in name.attributes){
                    args.pointEstimate[[iAttr]] <- attr(pointEstimate, iAttr)
                }
                
                outIID <- do.call(iidATE, args.pointEstimate)
            }else{ ## ignore that the nuisance parameters have been estimated
                outIID <- attr(pointEstimate, "iid.ate")
            }
            bootseeds <- NULL
            boot.object <- NULL

            ## display end
            if (verbose>1){ 
                cat("done\n")
            }

   
                                        # }}}
        }
    } else{
        outIID <- NULL
        boot.object <- NULL
        bootseeds <- NULL
                                        # }}}
    }
                                        # {{{ output object
    ## ** export
    if(!is.null(strata)){ ## when computing the average risk over strata
        setnames(pointEstimate$meanRisk,
                 old = "Treatment",
                 new = strata)
        setnames(pointEstimate$riskComparison,
                 old = c("Treatment.A","Treatment.B"),
                 new = paste0(strata,c(".A",".B")))
    }
    
    out <- list(meanRisk = pointEstimate$meanRisk,
                riskComparison = pointEstimate$riskComparison,
                iid = outIID,
                treatment = treatment,
                contrasts = contrasts,
                times = times,
                se = se,
                TD = TD,
                B = B,
                band = band,
                boot = boot.object,
                seeds = bootseeds)

  
    class(out) <- c("ate",class(object.event))
    if(se || band){
        if (verbose>1){ ## display
            cat(" - Confidence intervals: ")
        }
        ## FIXME: odd to have confint read and return everything 
        try.ci <- try(this <- stats::confint(out),silent=TRUE)
        if (class(try.ci)[1]=="try-error"){
        }else{
            out <- this
        }
        if(iid == FALSE){
            out$iid <- NULL
        }
        if (verbose>1){ ## display
            cat("done\n")
        }
    }

    return(out)

}

## * ate_initArgs
ate_initArgs <- function(object.event,
                         object.treatment,
                         object.censor,
                         landmark,
                         formula,
                         data,
                         cause,
                         times,
                         handler,
                         product.limit){

    ## ** user-defined arguments
    ## times
    if(inherits(object.event,"glm") && missing(times)){
        times <- NA
    }
    ## handler
    handler <- match.arg(handler, c("foreach","mclapply","snow","multicore"))

    ## ** fit regression model when user specifies formula and extract formula
    if(inherits(object.event,"formula")){
        myformula.event <- object.event

        if(grep("Hist",object.event)){
            object.event <- CSC(myformula.event, data = data, surv.type = "survival")
        }else{
            object.event <- glm(myformula.event, data = data, family = stats::binomial(link = "logit"))
        }
    }else{
        if(inherits(object.event,"glm") || inherits(object.event,"CauseSpecificCox")){
            myformula.event <- stats::formula(object.event)
        }else if(inherits(object.event,"coxph") || inherits(object.event,"cph") ||inherits(object.event,"phreg")){
            myformula.event <- coxFormula(object.event)
        }

    }

    if(!missing(object.treatment) && inherits(object.treatment,"formula")){
        myformula.treatment <- object.treatment
        object.treatment <- glm(myformula.treatment, data = data, family = stats::binomial(link = "logit"))        
    }else if(!missing(object.treatment) && inherits(object.treatment,"glm")){
        myformula.treatment <- stats::formula(object.treatment)
    }

    if(inherits(object.censor,"formula")){
        myformula.censor <- object.censor
        object.censor <- rms::cph(myformula.censor, data = data, x = TRUE, y = TRUE)
    }else if(inherits(object.censor,"coxph") || inherits(object.censor,"cph") || inherits(object.censor,"phreg")){
        myformula.censor <- coxFormula(object.censor)
    }

    
    ## ** deduced from user defined arguments

    ## censoring
    if(inherits(object.censor,"coxph") || inherits(object.censor,"cph") || inherits(object.censor,"phreg")){
        censoringMF <- coxModelFrame(object.censor)
        n.censor <- sapply(times, function(t){sum((censoringMF$status == 1) * (censoringMF$stop <= t))})
    }else{
        n.censor <- 0
    }
    
    ## level.censoring
    if(!is.null(object.censor) && (inherits(object.censor,"coxph") || inherits(object.censor,"cph") || inherits(object.censor,"phreg"))){
        censorVar <- SurvResponseVar(myformula.censor)
        censorVar.status <- censorVar$status
        censorVar.time <- censorVar$time
        if(censorVar.status %in% names(data)){
            iCandidates <- data[[censorVar.status]][coxModelFrame(object.censor)$status==1]
            level.censoring <- unique(iCandidates)[1]            
        }else{
            level.censoring <- NA
        }
    }else{
        level.censoring <- NA
        censorVar.status <- NA
        censorVar.time <- NA
    }

    ## event
    if(inherits(object.event,"glm")){
        eventVar.time <- as.character(NA)
        eventVar.status <- all.vars(myformula.event)[1]
        type.multistate <- "NA"
    }else if(inherits(object.event,"coxph") || inherits(object.event,"cph") || inherits(object.event,"phreg")){
        responseVar <- SurvResponseVar(myformula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        type.multistate <- "survival"
    }else if(inherits(object.event,"CauseSpecificCox")){
        causeNNA <- na.omit(c(cause,1))[1]
        responseVar <- SurvResponseVar(myformula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        type.multistate <- "competing.risks"
    }else{
        eventVar.time <- object.event[1]
        eventVar.status <- object.event[2]
        object.event <- NULL
        
        all.states_nC <- setdiff(unique(data[[eventVar.status]]),level.censoring)
        type.multistate <- if(length(all.states_nC)>1){"competing.risks"}else{"survival"}
    }

    ## treatment
    if(missing(object.treatment)){
        treatment <- NULL
        object.treatment <- NULL
        object.censor <- NULL
    }else  if(inherits(object.treatment,"glm")){
        treatment <- all.vars(myformula.treatment)[1]
    }else{
        treatment <- object.treatment
        object.treatment <- NULL
        object.censor <- NULL
    }

    ## cause
    if(is.na(cause) && (eventVar.status %in% names(data))){

        if(!is.null(object.event) && !inherits(object.event,"CauseSpecificCox")){
            event.call <- attr(SurvResponseVar(myformula.event)$status,"call")
        }else{
            event.call <- NULL
        }
            
        if(!is.null(event.call) && grepl("==",event.call)){
            event.call.rhs <- trimws(gsub(")","",strsplit(event.call,"==")[[1]][2]), which = "both")
            cause <- eval(parse(text = event.call.rhs))
        }else if(is.factor(data[[eventVar.status]]) && length(levels(data[[eventVar.status]]))==2){
            cause <- levels(data[[eventVar.status]])[2]
        }else if(all(data[[eventVar.status]] %in% c(0,1)) ){
            cause <- 1
        }
    }

    ## presence of time dependent covariates
    TD <- switch(class(object.event)[[1]],
                 "coxph"=(attr(object.event$y,"type")=="counting"),
                 "CauseSpecificCox"=(attr(object.event$models[[1]]$y,"type")=="counting"),
                 FALSE)
    if(TD){
        fct.pointEstimate <- ATE_TD
    }else{
        landmark <- NULL
        formula <- NULL
        fct.pointEstimate <- ATE_TI
    }

    n.obs <- NROW(data)
        
    ## estimator
    if(!is.null(object.event) && is.null(object.treatment)){
        if(TD){
            estimator <- "GformulaTD"
        }else{
            estimator <- "Gformula"
        }
    }else if(is.null(object.event) && !is.null(object.treatment)){
        if(any(n.censor>0)){
            estimator <- "IPTW,IPCW"
        }else{
            estimator <- "IPTW"
        }
    }else if(!is.null(object.event) && !is.null(object.treatment)){
        if(any(n.censor>0)){
            estimator <- "AIPTW,AIPCW"
        }else{
            estimator <- "AIPTW"
        }
    }else{
        estimator <- NA
    }

    #### which terms are used by the estimator
    ## term 1: F_1(\tau|A=a,W) or F_1(\tau|A=a,W) (1-1(A=a)/Prob[A=a|W])
    if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
        attr(estimator,"Gformula") <- TRUE
    }else{
        attr(estimator,"Gformula") <- FALSE
    }

    ## term 2 (treatment): 1(A=a)Y(tau)/Prob[A=a|W] or 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
        attr(estimator,"IPTW") <- TRUE
    }else{
        attr(estimator,"IPTW") <- FALSE
    }

    ## term 2 (censoring): 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW")){
        attr(estimator,"IPCW") <- TRUE
    }else{
        attr(estimator,"IPCW") <- FALSE
    }

    ## term 3: 1(A=a)Y(tau) \int () dM
    if(estimator %in% c("AIPTW,AIPCW")){
        attr(estimator,"integral") <- TRUE
    }else{
        attr(estimator,"integral") <- FALSE
    }

    ## ** product.limit
    if(is.null(product.limit)){
        product.limit <- switch(type.multistate,
                                "survival" = FALSE,
                                "competing.risks" = TRUE,
                                "NA" = NA)
    }
    
    ## ** output
    return(list(object.event = object.event,
                object.treatment = object.treatment,
                object.censor = object.censor,
                times = times,
                handler = handler,
                estimator = estimator,
                formula = formula,
                landmark = landmark,
                TD = TD,
                fct.pointEstimate = fct.pointEstimate,
                n.obs = n.obs,
                n.censor = n.censor,
                treatment = treatment,
                cause = cause,
                level.censoring = level.censoring,
                eventVar.time = eventVar.time,
                eventVar.status = eventVar.status,
                censorVar.time = censorVar.time,
                censorVar.status = censorVar.status,
                product.limit = product.limit
                ))
}

## * ate_checkArgs
ate_checkArgs <- function(object.event,
                          object.treatment,
                          object.censor,
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
                          TD,
                          n.censor,
                          level.censoring,
                          estimator,
                          eventVar.time,
                          eventVar.status,
                          censorVar.time,
                          censorVar.status){

    options <- riskRegression.options()
    
    ## ** times
    if(inherits(object.event,"glm") && length(times)!=1){
        stop("Argument \'times\' has no effect when using a glm object \n",
             "It should be set to NA \n")        
    }

    ## ** object.event
    if(!is.null(object.event)){
        candidateMethods <- paste("predictRisk",class(object.event),sep=".")
        if (all(match(candidateMethods,options$method.predictRisk,nomatch=0)==0)){
            stop(paste("Could not find predictRisk S3-method for ",class(object.event),collapse=" ,"),sep="")
        }

        if(estimator %in% c("AIPTW","AIPTW,AIPCW") && inherits(object.event,"glm")){
            warnings("It is unclear whether the current implementation of the double robust estimator is valid for logistic models.\n")
        }
    }
    
    ## ** object.treatment    
    if(!is.null(object.treatment)){
        if(!inherits(object.treatment,"glm") || object.treatment$family$family!="binomial"){
            stop("Argument \'treatment\' must be a logistic regression\n",
                 " or a character variable giving the name of the treatment variable. \n")
        }
    }

    ## ** object.censor
    if(estimator %in% c("AIPTW,AIPCW","IPTW,IPCW")){
        ## require censoring model
        if(is.null(object.censor)){
            stop("Argument \'censor\' must not be NULL for ",estimator," estimators in presence of censoring")
        }
        if(!inherits(object.censor,"coxph") && !inherits(object.censor,"cph") && !inherits(object.censor,"phreg")){
            stop("Argument \'censor\' must be a Cox model \n")
        }
        if(!identical(censorVar.status,eventVar.status)){
            if(sum(diag(table(data[[censorVar.status]],data[[eventVar.status]])))!=NROW(data)){
                stop("The status variables in object.event and object.censor are inconsistent \n")
            }
        }
        causeANDcensor <- intersect(which(data[[censorVar.status]]==level.censoring),
                                    which(data[[eventVar.status]]==cause))
        if(length(causeANDcensor)>0){
            stop("The status variables in object.event and object.censor are inconsistent \n",
                 "Some observations are both censored and have the event of interest \n")        
        }
        
        if(!identical(censorVar.time,eventVar.time)){
            stop("The time variable should be the same in object.event and object.censor \n")
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
        if(!is.null(object.event) && is.null(object.event$call)){
            stop("Argument \'event\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        if(!is.null(object.treatment) && is.null(object.treatment$call)){
            stop("Argument \'treatment\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        if(!is.null(object.censor) && is.null(object.censor$call)){
            stop("Argument \'censor\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
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

    ## ** functional delta method
    if(B==0 && (se|band|iid)){
        if(!is.null(landmark)){
            stop("Calculation of the standard errors via the functional delta method not implemented for time dependent covariates \n")
        }

        if(!is.null(object.event)){
            candidateMethods <- paste("predictRiskIID",class(object.event),sep=".")
            if (all(match(candidateMethods,options$method.predictRiskIID,nomatch=0)==0)){
                stop(paste("Could not find predictRiskIID S3-method for ",class(object.event),collapse=" ,"),"\n",
                     "Functional delta method not implemented for this type of object \n",
                     "Set argument \'B\' to a positive integer to use a boostrap instead \n",sep="")
            }

            if(coxN(object.event)[1]!=NROW(data)){
                stop("Argument \'event\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                     "(discrepancy found in number of rows) \n")
            }
        }

        if(!is.null(object.treatment) && coxN(object.treatment)!=NROW(data)){ ## note: only check that the datasets have the same size
            stop("Argument \'treatment\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                 "(discrepancy found in number of rows) \n")
        }

        if(!is.null(object.censor) && coxN(object.censor)!=NROW(data)){ ## note: only check that the datasets have the same size
                stop("Argument \'censor\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                     "(discrepancy found in number of rows) \n")
        }

    }

    
    ## ** treatment
    if(!is.null(treatment)){        
        if(treatment %in% names(data) == FALSE){
            stop("The data set does not seem to have a variable ",treatment," (argument: object.treatment). \n")
        }
        if(is.numeric(data[[treatment]])){
            stop("The treatment variable must be a factor variable. \n",
                 "Convert treatment to factor, re-fit the object using this new variable and then call ate. \n")
        }
    }else if(is.null(strata)){
        stop("The treatment variable must be specified using the argument \'treatment\' \n")
    }

    ## ** strata
    if(!is.null(strata)){
        if(estimator != "Gformula"){
            stop("Argument \'strata\' only compatible with the G-formula estimator \n")
        }
        if(any(strata %in% names(data) == FALSE)){
            stop("The data set does not seem to have a variable \"",paste0(strata, collapse = "\" \""),"\" (argument: strata). \n")
        }
        if(length(strata) != 1){
            stop("Argument strata should have length 1. \n")
        }
        if(TD == TRUE){
            stop("Landmark analysis is not available when argument strata is specified. \n")
        }
    }

    ## ** status
    if(eventVar.status %in% names(data) == FALSE){
        stop("The data set does not seem to have a variable ",eventVar.status," (argument: object.event[2]). \n")
    }
    if(is.na(cause)){
        stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
             "Please specify the argument \'cause\'. \n")
    }
    if(cause %in% data[[eventVar.status]] == FALSE){
        stop("Value of the argument \'cause\' not found in column \"",eventVar.status,"\" \n")
    }
    
    ## ** event time
    if(!is.na(eventVar.time)){
        if(eventVar.time %in% names(data) == FALSE){
            stop("The data set does not seem to have a variable ",eventVar.time," (argument: object.event[1]). \n")
        }
        if(is.na(cause)){
            stop("Argument \'cause\' not specified\n")
        }

        if(!is.null(treatment)){
            freq.event <- tapply(data[[eventVar.status]], data[[treatment]], function(x){mean(x==cause)})
            count.event <- tapply(data[[eventVar.status]], data[[treatment]], function(x){sum(x==cause)})
        }else{
            freq.event <- tapply(data[[eventVar.status]], data[[strata]], function(x){mean(x==cause)})
            count.event <- tapply(data[[eventVar.status]], data[[strata]], function(x){sum(x==cause)})
        }
        if(any(freq.event < 0.01) || any(count.event < 5)  ){
            warning("Rare event, possible violation of the positivity assumption \n")
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
    
    ## ** output
    return(TRUE)
}


            

#----------------------------------------------------------------------
### ate.R ends here
