### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: sep 23 2020 (16:21) 
##           By: Brice Ozenne
##     Update #: 2061
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * ate (documentation)
#' @title Average Treatment Effects Computation 
#' @description Use the g-formula/IPTW/double robust estimator to estimate the average treatment
#'     effect based on Cox regression with or without competing risks.
#' @name ate
#' 
#' @param event Outcome model which describes how the probability of experiencing a terminal event depends
#'     on treatment and covariates. The object carry its own call and
#'     have a \code{predictRisk} method. See examples.
#' @param treatment Treatment model which describes how the probability of being allocated to a treatment group depends
#'     on covariates. The object must be a \code{glm} object (logistic regression) or the name of the treatment variable.
#'     See examples.
#' @param censor Censoring model which describes how the probability of being censored depends
#'     on treatment and covariates. The object must be a \code{coxph} or \code{cph} object. See examples.
#' @param data [data.frame or data.table] Data set in which to evaluate risk predictions
#' based on the outcome model
#' @param data.index [numeric vector] Position of the observation in argument data relative to the dataset used to obtain the argument event, treatment, censor.
#' Only necessary for the standard errors when computing the Average Treatment Effects on a subset of the data set. 
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
#' @param estimator [character] The type of estimator used to compute the average treatment effect. 
#' Can be \code{"G-formula"}, \code{"IPTW"}, or \code{"AIPTW"}.
#' When using \code{estimator="G-formula"}, a model for the outcome should be provided (argument event).
#' When using \code{estimator="IPTW"}, a model for the treatment should be provided (argument treatment), as well as for the censoring (if any, argument censor).
#' When using \code{estimator="AIPTW"} (double robust estimator), a model for the outcome and the treatment should be provided (argument event and treatment), as well as for the censoring (if any, argument censor).
#' @param known.nuisance [logical] If \code{FALSE} the uncertainty related to the estimation of the nuisance parameters is ignored.
#' This greatly simplifies computations but requires to use a double robust estimator.
#' The resulting standard error is known to be consistent when all event, treatment, and censoring models are valid.
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
#' @param ... passed to predictRisk
#'
#' @author Brice Ozenne \email{broz@@sund.ku.dk}
#' and Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#'
#' @seealso
#' \code{\link{confint.ate}} to compute confidence intervals/bands.
#' \code{\link{autoplot.ate}} to display the average risk.

## * ate (examples)
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
#' ## standard error computed using the influence function
#' ## confidence intervals / p-values based on asymptotic results
#' ateFit1a <- ate(fit, data = dtS, treatment = "X1", times = 5:8)
#' summary(ateFit1a)
#' 
#' \dontrun{
#' ## same as before with in addition the confidence bands / adjusted p-values
#' ## (argument band = TRUE)
#' ateFit1b <- ate(fit, data = dtS, treatment = "X1", times = 5:8,
#'                 band = TRUE)
#' summary(ateFit1b)
#'
#' ## by default bands/adjuste p-values computed separately for each treatment modality
#' print(confint(ateFit1b, band = 1), quantile = TRUE)
#' ## adjustment over all treatment and times using the band argument of confint
#' print(confint(ateFit1b, band = 2), quantile = TRUE)
#' 
#' ## confidence intervals / p-values computed using 1000 boostrap samples
#' ## (argument se = TRUE and B = 1000) 
#' ateFit1c <- ate(fit, data = dtS, treatment = "X1",
#'                 times = 5:8, se = TRUE, B = 1000)
#' ## NOTE: for real applications 1000 bootstrap samples is not enough 
#'
#' ## same but using 2 cpus for generating and analyzing the boostrap samples
#' ## (parallel computation, argument mc.cores = 2) 
#' ateFit1d <- ate(fit, data = dtS, treatment = "X1",
#'                 times = 5:8, se = TRUE, B = 1000, mc.cores = 2)
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
#' ## with confidence intervals
#' ateFit1b <- ate(fit, data = dtB, treatment = "X1",
#'                times = 5) ## just for having a nice output not used in computations
#' ateFit1b
#' 
#' ## using the lava package
#' ateLava <- estimate(fit, function(p, data){
#' a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X2"] ;
#' R.X11 <- expit(a + b + c * data[["X2"]])
#' R.X10 <- expit(a + c * data[["X2"]])
#' list(risk0=R.X10,risk1=R.X11,riskdiff=R.X11-R.X10)},
#' average=TRUE)
#' ateLava
#'
#' ateFit1b$meanRisk
#' }
#' 
#' #### Competing risks settings               ####
#' #### ATE with cause specific Cox regression ####
#'
#' ## generate data
#' n <- 500
#' set.seed(10)
#' dt <- sampleData(n, outcome="competing.risks")
#' dt$X1 <- factor(rbinom(n, prob = c(0.2,0.3) , size = 2), labels = paste0("T",0:2))
#'
#' ## estimate cause specific Cox model
#' fitCR <-  CSC(Hist(time,event)~ X1+X8,data=dt,cause=1)
#' 
#' ## compute the ATE at times 1, 5, 10 using X1 as the treatment variable
#' ateFit2a <- ate(fitCR, data = dt, treatment = "X1", times = c(1,5,10),
#'                 cause = 1, se = TRUE, band = TRUE)
#' ateFit2a
#' 
#' #### Double robust estimator ####
#' \dontrun{
#' ## generate data
#' n <- 500
#' set.seed(10)
#' dt <- sampleData(n, outcome="competing.risks")
#' dt$time <- round(dt$time,1)
#' dt$X1 <- factor(rbinom(n, prob = c(0.4) , size = 1), labels = paste0("T",0:1))
#'
#' ## working models
#' m.event <-  CSC(Hist(time,event)~ X1+X2+X3+X5+X8,data=dt)
#' m.censor <-  coxph(Surv(time,event==0)~ X1+X2+X3+X5+X8,data=dt, x = TRUE, y = TRUE)
#' m.treatment <-  glm(X1~X2+X3+X5+X8,data=dt,family=binomial(link="logit"))
#'
#' ## prediction + average
#' ateRobust <- ate(event = m.event,
#'                  treatment = m.treatment,
#'                  censor = m.censor,
#'                  data = dt, times = 5:10, 
#'                  cause = 1, band = TRUE)
#' 
#' ## compare various estimators
#' ateRobust3 <- ate(event = m.event,
#'                  treatment = m.treatment,
#'                  censor = m.censor,
#'                  estimator = c("GFORMULA","IPTW","AIPTW"),
#'                  data = dt, times = c(5:10), 
#'                  cause = 1, se = TRUE)
#' ateRobust3$meanRisk[,.(meanRisk.Gformula,meanRisk.IPTW,meanRisk.AIPTW)]
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
#'         landmark = c(0,30,60,90), cause = 1, B = 10, se = 1,
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
#' cox1TD <- coxph(Surv(start,time, status,type="counting") ~ X3+X5+X6+X8,
#'                 data=d, x = TRUE)
#' resTD1 <- ate(cox1TD,formula=Hist(entry=start,time=time,event=status)~1,
#'         data = d, treatment = "X3", contrasts = NULL,
#'         times=.5,verbose=1,
#'         landmark = c(0,0.5,1), B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resTD1
#' ## account for competing risks
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
                data.index = NULL,
                formula,
                estimator = NULL,
                strata = NULL,
                contrasts = NULL,
                times,
                cause = NA,
                landmark,
                se = TRUE,
                iid = (B == 0) && (se || band),
                known.nuisance = FALSE,
                band = FALSE,
                B = 0,
                seed,
                handler = "foreach",
                mc.cores = 1,
                cl = NULL,
                verbose = TRUE,
                ...){

    dots <- list(...)

    ## ** initialize arguments
    if(is.data.table(data)){
        data <- data.table::copy(data)
    }else{
        data <- data.table::as.data.table(data)
    }

    init <- ate_initArgs(object.event = event,
                         object.treatment = treatment,
                         object.censor = censor,
                         formula = formula,
                         landmark = landmark,
                         mydata = data,
                         data.index = data.index,
                         estimator = estimator,
                         times = times,
                         cause = cause,
                         handler = handler,
                         product.limit = dots$product.limit,
                         store.iid = dots$store.iid)

    object.event <- init$object.event
    object.treatment <- init$object.treatment
    object.censor <- init$object.censor
    data.index <- init$data.index
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
    level.states <- init$level.states
    censorVar.time <- init$censorVar.time
    censorVar.status <- init$censorVar.status
    dots$product.limit <- init$product.limit
    store.iid <- init$store.iid

    return.iid.nuisance <- (se || band || iid) && (B==0) && (known.nuisance == FALSE)
    ## method to compute the iid for the ate:
    ## - 1) implicit average: call predictRisk with average.iid=TRUE which performs the average
    ## - 2) explicit average: call predictRisk with iid=TRUE and then perform the average (slower)
    if("method.iid" %in% names(dots)){ 
        method.iid <- dots$method.iid
        dots$method.iid <- NULL
    }else{
        method.iid <- 1 ## implicit
    }
    
    ## ** check consistency of the user arguments
    check <- ate_checkArgs(object.event = object.event,
                           object.censor = object.censor,
                           object.treatment = object.treatment,
                           mydata = data,
                           formula = formula,
                           treatment = treatment,
                           strata = strata,
                           contrasts = contrasts,
                           times = times,
                           cause = cause,
                           landmark = landmark,
                           se = se,
                           iid = iid,
                           data.index = data.index,
                           return.iid.nuisance = return.iid.nuisance,
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
                           level.states = level.states,
                           estimator = estimator,
                           eventVar.time = eventVar.time,
                           eventVar.status = eventVar.status,
                           censorVar.time = censorVar.time,
                           censorVar.status = censorVar.status,
                           n.obs = n.obs
                           )

    ## ** Prepare
    testE.Cox <- inherits(object.event,"coxph") || inherits(object.event,"cph") || inherits(object.event,"phreg")
    testE.CSC <- inherits(object.event,"CauseSpecificCox")
    testE.glm <- inherits(object.event,"glm")

    if(!is.null(strata)){
        if(!is.factor(data[[strata]])){
            data[[strata]] <- factor(data[[strata]])
        }   
        if(is.null(contrasts)){
            contrasts <- levels(data[[strata]])
        }
        levels <- levels(data[[strata]]) ## necessary to get the correct factor format in case that not all levels are in contrast
    }else{
        if(!is.factor(data[[treatment]])){
            data[[treatment]] <- factor(data[[treatment]])
        }   
        if(is.null(contrasts)){
            contrasts <- levels(data[[treatment]])
            ## if (length(contrasts)>50) warning("Treatment variable has more than 50 levels.\nIf this is not a mistake,
            ## you should use the argument `contrasts'.")
        }
        levels <- levels(data[[treatment]]) ## necessary to get the correct factor format in case that not all levels are in contrast
    }

    ## object to be exported
    out <- list(meanRisk = NULL,
                diffRisk = NULL,
                ratioRisk = NULL,
                iid = NULL,
                boot = NULL,
                estimator = estimator,
                eval.times = times, 
                n = n.obs, 
                variables = c(time = eventVar.time, event = eventVar.status,
                              treatment = if(is.null(treatment)){NA}else{treatment},
                              strata = if(is.null(strata)){NA}else{strata}),
                theCause = cause,
                contrasts = contrasts,
                computation.time = list("point" = NA, "iid" = NA, "bootstrap" = NA), 
                inference = data.frame("se" = se, "iid" = iid, "band" = band, "B" = B, "ci" = FALSE, "p.value" = FALSE,
                                       "conf.level" = NA, "alternative" = NA,
                                       "bootci.method" = NA,
                                       "method.band" = NA, "n.sim" = NA)
                )
    attr(out$theCause,"cause") <- level.states
    attr(out$theCause,"level.censoring") <- level.censoring
    attr(out$variables,"range.time") <- range(data[[eventVar.time]])
    if(!is.null(strata)){
        attr(out$eval.times,"n.at.risk") <- dcast(cbind(times, data[, .(pc = sapply(times, function(t){sum(.SD[[eventVar.time]]>=t)})), by = strata]),
                                                  value.var = "pc", formula = as.formula(paste0(strata,"~times")))
    }else{
        attr(out$eval.times,"n.at.risk") <- dcast(cbind(times, data[, .(pc = sapply(times, function(t){sum(.SD[[eventVar.time]]>=t)})), by = treatment]),
                                                  value.var = "pc", formula = as.formula(paste0(treatment,"~times")))
    }
    attr(out$eval.times,"n.censored") <- n.censor
    class(out) <- c("ate")

    ## ** Overview of the calculations
    if(verbose>1/2){
        print(out)
        cat("\n Processing\n")
    }

    ## ** Compute influence function relative to each Cox model (if not already done)
    if(return.iid.nuisance){
        if(verbose>1/2){
            cat(" - Prepare influence function:")
        }

        if(attr(estimator,"GFORMULA") && (testE.Cox || testE.CSC) && identical(is.iidCox(object.event),FALSE)){
            if(verbose>1/2){cat(" outcome")}
            object.event <- iidCox(object.event, tau.max = max(times), store.iid = store.iid, return.object = TRUE)
        }

        if(attr(estimator,"IPCW")){
            if(identical(is.iidCox(object.censor),FALSE)){
                if(verbose>1/2){cat(" censoring")}
                object.censor <- iidCox(object.censor, store.iid = store.iid, tau.max = max(times))
            }
        }
        if(verbose>1/2){cat(" done \n")}
    }
    
    ## ** Point estimate
    if(verbose>1/2){
        cat(" - Point estimation:")
    }
    args.pointEstimate <- c(list(object.event = object.event,
                                 object.censor = object.censor,
                                 object.treatment = object.treatment,
                                 mydata=data,
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
                                 return.iid.nuisance = return.iid.nuisance,
                                 data.index = data.index,
                                 method.iid = method.iid),
                            dots)
    if (TD){       
        args.pointEstimate <- c(args.pointEstimate,list(formula=formula))
    }
    ## note: system.time() seems to slow down the execution of the function, this is why Sys.time is used instead.
    tps1 <- Sys.time()

    pointEstimate <- do.call(fct.pointEstimate, args.pointEstimate)
    
    tps2 <- Sys.time()

    out$meanRisk <- pointEstimate$meanRisk
    out$diffRisk <- pointEstimate$diffRisk
    out$ratioRisk <- pointEstimate$ratioRisk
    out$computation.time[["point"]] <- tps2-tps1

    if(verbose>1/2){cat(" done \n")}

    ## ** Confidence intervals
    if(se || band || iid){
        if (TD){
            key1 <- c("treatment","landmark")
            key2 <- c("treatment.A","treatment.B","landmark")
        }
        else{
            key1 <- c("treatment","time")
            key2 <- c("treatment.A","treatment.B","time")
        }

        if(B>0){
                                        # {{{ Bootstrap
            ## *** Bootstrap
            ## display start
            if (verbose>1/2){ 
                cat(" - Non-parametric bootstrap using ",B," samples and ",mc.cores," core",if(mc.cores>1){"s"},"\n", sep ="")
                cat("                            (expected time: ",round(as.numeric(out$computation.time$point)*B/mc.cores,2)," seconds)\n", sep = "")
            }

            name.estimate <- c(paste("mean",out$meanRisk$estimator,out$meanRisk$time,out$meanRisk$landmark,out$meanRisk$treatment,sep="."),
                               paste("difference",out$diffRisk$estimator,out$diffRisk$time,out$diffRisk$landmark,out$diffRisk$A,out$diffRisk$B,sep="."),
                               paste("ratio",out$ratioRisk$estimator,out$ratioRisk$time,out$ratioRisk$landmark,out$ratioRisk$A,out$ratioRisk$B,sep="."))
            estimate <- setNames(c(out$meanRisk$estimate,out$diffRisk$estimate,out$ratioRisk$estimate), name.estimate)#

            ## run
            resBoot <- calcBootATE(args = args.pointEstimate,
                                   n.obs = n.obs["data"],
                                   fct.pointEstimate = fct.pointEstimate,
                                   name.estimate = name.estimate,
                                   handler = handler,
                                   B = B,
                                   seed = seed,
                                   mc.cores = mc.cores,
                                   cl = cl,
                                   verbose = verbose)

            ## store
            out$boot <- list(t0 = estimate,
                             t = resBoot$boot,
                             R = B,
                             data = data,
                             seed = resBoot$bootseeds,
                             statistic = NULL,
                             sim = "ordinary",
                             call = quote(boot(data = XX, statistic = XX, R = XX)),
                             stype = "i",
                             strata = rep(1,n.obs["data"]),
                             weights = rep(1/n.obs["data"],n.obs["data"]),
                             pred.i = NULL,  ## Omitted if m is 0 or sim is not "ordinary" (from doc of boot::boot)
                             L = NULL, ## only used when sim is "antithetic" (from doc of boot::boot)
                             ran.gen = NULL, ## only used when sim is "parametric" (from doc of boot::boot)
                             mle = NULL ## only used when sim is "parametric" (from doc of boot::boot)
                             )
            class(out$boot) <- "boot"
            tps3 <- Sys.time()
            out$computation.time[["bootstrap"]] <- tps3-tps2

            ## display end
            if (verbose>1/2){ 
                cat("\n")
            }

                                        # }}}
        } else {
                                        # {{{ compute standard error and quantiles via the influence function
            ## *** Delta method
            ## display start
            if (verbose>1/2){ 
                cat(" - Decomposition iid: ")
            }

            if(return.iid.nuisance && (attr(estimator,"IPTW") == TRUE)){
                ## compute iid decomposition relative to the nuisance parameters

                ## add pre-computed quantities
                args.pointEstimate[names(pointEstimate$store)] <- pointEstimate$store

                ## compute extra term and assemble
                if(identical(method.iid,2)){
                    out$iid <- do.call(iidATE2, args.pointEstimate)
                }else{
                    out$iid <- do.call(iidATE, args.pointEstimate)
                }
            }else{
                ## otherwise no extra computation required
                out$iid <- list(GFORMULA = pointEstimate$store$iid.GFORMULA,
                                IPTW = pointEstimate$iid.IPTW,
                                AIPTW = pointEstimate$iid.AIPTW)
            }
            tps3 <- Sys.time()
            out$computation.time[["iid"]] <- tps3-tps2

            ## display end
            if (verbose>1/2){ 
                cat("done\n")
            }

   
                                        # }}}
        }
    }
                                        # {{{ output object
    ## ** statistical inference
    if(se || band){
        if (verbose>1/2){ ## display
            if(se && band){
                cat(" - Confidence intervals / bands: ")
            }else if(se){
                cat(" - Confidence intervals: ")
            }else if(band){
                cat(" - Confidence bands: ")
            }
        }
        ## FIXME: odd to have confint read and return everything
        out <- stats::confint(out)
        if(iid == FALSE){
            out$iid <- NULL
        }
        if (verbose>1/2){ ## display
            cat("done\n")
        }
    }

    ## ** export
    return(out)

}

## * ate_initArgs
ate_initArgs <- function(object.event,
                         object.treatment,
                         object.censor,
                         landmark,
                         formula,
                         estimator,
                         mydata, ## instead of data to avoid confusion between the function data and the dataset when running the bootstrap
                         data.index,
                         cause,
                         times,
                         handler,
                         product.limit,
                         store.iid){

    ## ** user-defined arguments
    ## times
    if(inherits(object.event,"glm") && missing(times)){
        times <- NA
    }
    ## handler
    handler <- match.arg(handler, c("foreach","mclapply","snow","multicore"))
    ## data.index
    if(is.null(data.index)){data.index <- 1:NROW(mydata)}
    
    ## ** fit a regression model when the user specifies a formula
    ## Event
    if(missing(object.event)){
        formula.event <- NULL
        model.event <- NULL
    }else if(inherits(object.event,"formula")){ ## formula
        formula.event <- object.event
        if(any(grepl("Hist(",object.event, fixed = TRUE))){
            model.event <- do.call(CSC, args = list(formula = object.event, data = mydata))
        }else if(any(grepl("Surv(",object.event, fixed = TRUE))){
            model.event <- do.call(rms::cph, args = list(formula = object.event, data = mydata, x = TRUE, y = TRUE))
        }else{
            model.event <- do.call(stats::glm, args = list(formula = object.event, data = mydata, family = stats::binomial(link = "logit")))
        }
    }else if(is.list(object.event) && all(sapply(object.event, function(iE){inherits(iE,"formula")}))){   ## list of formula
        formula.event <- object.event
        model.event <- CSC(object.event, data = mydata)
    }else if(inherits(object.event,"glm") ||inherits(object.event,"wglm") || inherits(object.event,"CauseSpecificCox")){ ## glm / CSC
        formula.event <- stats::formula(object.event)
        model.event <- object.event
    }else if(inherits(object.event,"coxph") || inherits(object.event,"cph") ||inherits(object.event,"phreg")){ ## Cox
        formula.event <- coxFormula(object.event)
        model.event <- object.event
    }else{ ## character (time,event) i.e. no model
        formula.event <- NULL
        model.event <- NULL
    }

    ## Treatment
    if(missing(object.treatment)){
        formula.treatment <- NULL
        model.treatment <- NULL
    }else if(inherits(object.treatment,"formula")){
        model.treatment <- do.call(stats::glm, args = list(formula = object.treatment, data = mydata, family = stats::binomial(link = "logit")))        
        formula.treatment <- object.treatment
    }else if(inherits(object.treatment,"glm")){
        formula.treatment <- stats::formula(object.treatment)
        model.treatment <- object.treatment
    }else{
        formula.treatment <- NULL
        model.treatment <- NULL
    }

    ## Censoring
    if(missing(object.censor)){
        formula.censor <- NULL
        model.censor <- NULL
    }else if(inherits(object.censor,"formula")){
        formula.censor <- object.censor
        model.censor <- do.call(rms::cph, args = list(formula = object.censor, data = mydata, x = TRUE, y = TRUE))
    }else if(inherits(object.censor,"coxph") || inherits(object.censor,"cph") || inherits(object.censor,"phreg")){
        formula.censor <- coxFormula(object.censor)
        model.censor <- object.censor
    }else{
        formula.censor <- NULL
        model.censor <- NULL
    }

    ## ** identify event, treatment and censoring variables
    if(inherits(model.censor,"coxph") || inherits(model.censor,"cph") || inherits(model.censor,"phreg")){
        censoringMF <- coxModelFrame(model.censor)
        test.censor <- censoringMF$status == 1
        n.censor <- sapply(times, function(t){sum(test.censor * (censoringMF$stop <= t))})

        info.censor <- SurvResponseVar(coxFormula(model.censor))
        censorVar.status <- info.censor$status
        censorVar.time <- info.censor$time
        level.censoring <- unique(mydata[[censorVar.status]][model.censor$y[data.index,2]==1])
    }else{ ## G-formula or IPTW (no censoring)
        censorVar.status <- NA
        censorVar.time <- NA
        
        
        if(inherits(model.event,"CauseSpecificCox")){
            test.censor <- model.event$response[,"status"] == 0        
            n.censor <- sapply(times, function(t){sum(test.censor * (model.event$response[,"time"] <= t))})
            level.censoring <- attr(model.event$response,"cens.code")
        }else if(inherits(model.event,"coxph") || inherits(model.event,"cph") || inherits(model.event,"phreg")){ 
            censoringMF <- coxModelFrame(model.event)
            test.censor <- censoringMF$status == 0
            n.censor <- sapply(times, function(t){sum(test.censor * (censoringMF$stop <= t))})
            level.censoring <- 0
        }else if(inherits(model.event,"wglm")){
            n.censor <- object.event$n.censor
            level.censoring <- setdiff(unique(object.event$data[[object.event$var.outcome]]), object.event$causes)
        }else{
            n.censor <- rep(0,length(times))
        }
        if(all(n.censor==0)){level.censoring <- NA}
    }

    ## treatment
    if(missing(object.treatment) || is.null(object.treatment)){
        treatment <- NULL
    }else if(!is.null(formula.treatment)){
        treatment <- all.vars(formula.treatment)[1]
    }else if(is.character(object.treatment)){
        treatment <- object.treatment
    }
    
    ## event
    if(inherits(model.event,"CauseSpecificCox")){
        responseVar <- SurvResponseVar(formula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        if(is.na(cause)){
            cause <- model.event$theCause
        }
        if(is.null(product.limit)){product.limit <- TRUE}
        level.states <- model.event$cause
    }else if(inherits(model.event,"coxph") || inherits(model.event,"cph") || inherits(model.event,"phreg")){
        responseVar <- SurvResponseVar(formula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        if(is.na(cause)){ ## handle Surv(time,event > 0) ~ ...
            if(any(model.event$y[data.index,NCOL(model.event$y)]==1)){
                cause <- unique(mydata[[eventVar.status]][model.event$y[data.index,NCOL(model.event$y)]==1])
            }
        }
        if(is.null(product.limit)){product.limit <- FALSE}
        level.states <- 1
    }else if(inherits(object.event,"glm")){
        eventVar.time <- as.character(NA)
        eventVar.status <- all.vars(formula.event)[1]
        if(is.na(cause)){ ## handle I(Y > 0) ~ ...
            cause <- unique(mydata[[eventVar.status]][model.event$y[data.index]==1])
        }else{
            cause <- NA
        }
        if(is.null(product.limit)){product.limit <- NA}
        level.states <- unique(mydata[[eventVar.status]])
    }else if(inherits(object.event,"wglm")){
        eventVar.time <- object.event$var.time
        eventVar.status <- object.event$var.outcome
        if(is.na(cause)){ ## handle I(Y > 0) ~ ...
            cause <- object.event$theCause
        }
        level.states <- object.event$causes
        if(is.null(product.limit)){product.limit <- (length(level.states)>1)}
    }else{ ## no outcome model
        if(identical(names(object.event),c("time","status"))){
            eventVar.time <- object.event[1]
            eventVar.status <- object.event[2]
        }else if(identical(names(object.event),c("status","time"))){
            eventVar.status <- object.event[1]
            eventVar.time <- object.event[2]
        }else{            
            eventVar.time <- object.event[1]
            eventVar.status <- object.event[2]
        }

        level.states <- setdiff(unique(mydata[[eventVar.status]]), level.censoring)
        if(is.na(cause)){
            cause <- sort(level.states)[1]
        }
        if(is.null(product.limit)){product.limit <- NA}
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

    ## estimator
    if(is.null(estimator)){
        if(!is.null(model.event) && is.null(model.treatment)){
            if(TD){
                estimator <- "GFORMULATD"
            }else{
                estimator <- "GFORMULA"
            }
        }else if(is.null(model.event) && !is.null(model.treatment)){
            if(any(n.censor>0)){
                estimator <- "IPTW,IPCW"
            }else{
                estimator <- "IPTW"
            }
        }else if(!is.null(model.event) && !is.null(model.treatment)){
            if(any(n.censor>0)){
                estimator <- "AIPTW,AIPCW"
            }else{
                estimator <- "AIPTW"
            }
        }else{
            estimator <- NA
        }
    }else{
        estimator <- toupper(estimator)
        if(any(estimator == "IPTW") && any(n.censor>0)){
            estimator[estimator == "IPTW"] <- "IPTW,IPCW"
        }
        if(any(estimator == "AIPTW") && any(n.censor>0)){
            estimator[estimator == "AIPTW"] <- "AIPTW,AIPCW"
        }
    }

    #### which terms are used by the estimator
    if(any(estimator %in% c("GFORMULATD"))){
        attr(estimator,"TD") <- TRUE
    }else{
        attr(estimator,"TD") <- FALSE
    }

    ## term 1: F_1(\tau|A=a,W) or F_1(\tau|A=a,W) (1-1(A=a)/Prob[A=a|W])
    if(any(estimator %in% c("GFORMULA","GFORMULATD","AIPTW","AIPTW,AIPCW"))){
        attr(estimator,"GFORMULA") <- TRUE
    }else{
        attr(estimator,"GFORMULA") <- FALSE
    }

    ## term 2 (treatment): 1(A=a)Y(tau)/Prob[A=a|W] or 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(any(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW"))){
        attr(estimator,"IPTW") <- TRUE
    }else{
        attr(estimator,"IPTW") <- FALSE
    }

    ## term 2 (censoring): 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(any(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW"))){
        attr(estimator,"IPCW") <- TRUE
    }else{
        attr(estimator,"IPCW") <- FALSE
    }

    ## term 3: 1(A=a)Y(tau) \int () dM
    if(any(estimator %in% c("AIPTW,AIPCW"))){
        attr(estimator,"integral") <- !inherits(object.event,"wglm")
    }else{
        attr(estimator,"integral") <- FALSE
    }

    if(any(estimator %in% c("AIPTW","AIPTW,AIPCW"))){
        attr(estimator,"augmented") <- TRUE
    }else{
        attr(estimator,"augmented") <- FALSE
    }

    ## which estimates should be output?
    if(any(estimator %in% c("GFORMULA","GFORMULATD"))){
        attr(estimator,"export.GFORMULA") <- TRUE
    }else{
        attr(estimator,"export.GFORMULA") <- FALSE
    }
    if(any(estimator %in% c("IPTW","IPTW,IPCW"))){
        attr(estimator,"export.IPTW") <- TRUE
    }else{
        attr(estimator,"export.IPTW") <- FALSE
    }
    if(any(estimator %in% c("AIPTW","AIPTW,AIPCW"))){
        attr(estimator,"export.AIPTW") <- TRUE
    }else{
        attr(estimator,"export.AIPTW") <- FALSE
    }

    ## ** sample size
    n.obs <- c(data = NROW(mydata),
               model.event = if(!is.null(model.event)){stats::nobs(model.event)}else{NA},
               model.treatment = if(!is.null(model.treatment)){stats::nobs(model.treatment)}else{NA},
               model.censor = if(!is.null(model.censor)){stats::nobs(model.censor)}else{NA}
               )

    ## ** store.iid
    if(is.null(store.iid)){
        store.iid <- "minimal"
    }

    ## ** output
    return(list(object.event = model.event,
                object.treatment = model.treatment,
                object.censor = model.censor,
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
                level.states = level.states,
                eventVar.time = eventVar.time,
                eventVar.status = eventVar.status,
                censorVar.time = censorVar.time,
                censorVar.status = censorVar.status,
                product.limit = product.limit,
                data.index = data.index,
                store.iid = store.iid
                ))
}

## * ate_checkArgs
ate_checkArgs <- function(object.event,
                          object.treatment,
                          object.censor,
                          mydata,
                          formula,
                          treatment,
                          strata,
                          contrasts,
                          times,
                          cause,
                          landmark,
                          se,
                          iid,
                          data.index,
                          return.iid.nuisance,
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
                          level.states,
                          estimator,
                          eventVar.time,
                          eventVar.status,
                          censorVar.time,
                          censorVar.status,
                          n.obs){

    options <- riskRegression.options()

    ## ** times
    if(inherits(object.event,"glm") && length(times)!=1){
        stop("Argument \'times\' has no effect when using a glm object \n",
             "It should be set to NA \n")        
    }

    ## ** status
    if(eventVar.status %in% names(mydata) == FALSE){
        stop("The data set does not seem to have a variable ",eventVar.status," (argument: event[2]). \n")
    }
    if(inherits(object.event,"glm") && is.na(cause)){
        stop("Argument \'cause\' should not be specified when using a glm model in argument \'event\'")
    }
    if(length(cause)>1 || is.na(cause)){
        stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
             "Please specify the argument \'cause\'. \n")
    }
    if(cause %in% mydata[[eventVar.status]] == FALSE){
        stop("Value of the argument \'cause\' not found in column \"",eventVar.status,"\" \n")
    }

    ## ** estimator
    if(any(is.na(estimator))){
        stop("No model for the outcome/treatment/censoring has been specified. Cannot estimate the average treatment effect. \n")
    }
    
    valid.estimator <- c("GFORMULATD","GFORMULA","IPTW,IPCW","IPTW","AIPTW,AIPCW","AIPTW")
    if(any(estimator %in% valid.estimator == FALSE)){
        stop("Incorrect value(s) for argument \'estimator\': \"",paste(estimator[estimator %in% valid.estimator == FALSE], collapse = "\" \""),"\"\n",
             "Valid values: \"G-formula\", \"G-formulaTD\", \"IPTW\", \"AIPTW\" \n")
    }
    if(attr(estimator,"GFORMULA")){
        if(is.null(object.event)){
            stop("Using a G-formula/AIPTW estimator requires to specify a model for the outcome (argument \'event\') \n")
        }
    }
    if(attr(estimator,"IPTW")){
        if(is.null(object.treatment)){
            stop("Using a ITPW/AIPTW estimator requires to specify a model for the treatment allocation (argument \'treatment\') \n")
        }
    }
    if(attr(estimator,"IPCW")){
        if(is.null(object.censor)){
            stop("Using a ITPW/AIPTW estimator requires to specify a model for the censoring mechanism (argument \'censor\') in presence of right-censoring \n")
        }
    }
    
    ## ** object.event
    if(!is.null(object.event)){
        candidateMethods <- paste("predictRisk",class(object.event),sep=".")
        if (all(match(candidateMethods,options$method.predictRisk,nomatch=0)==0)){
            stop(paste("Could not find predictRisk S3-method for ",class(object.event),collapse=" ,"),sep="")
        }
        if((inherits(object.event,"wglm") || inherits(object.event,"CauseSpecificCox")) && (cause %in% object.event$causes == FALSE)){
            stop("Argument \'cause\' does not match one of the available causes: ",paste(object.event$causes,collapse=" "),"\n")
        }
        if(inherits(object.event,"wglm") && (cause != object.event$theCause)){
            stop("Argument \'cause\' does not match the one of \'object.event\' \n",
                 "Consider re-fitting the model with appropriate cause.")
        }
        if(any(is.na(level.censoring)) && any(na.omit(level.censoring) == cause)){
            stop("The cause of interest must differ from the level indicating censoring (in the outcome model) \n")
        } 
        if(any(is.na(coef(object.event)))){
            stop("Cannot handle missing values in the model coefficients (event) \n")
        }
    }
    
    ## ** object.treatment    
    if(!is.null(object.treatment)){
        if(!inherits(object.treatment,"multinom") && (!inherits(object.treatment,"glm") || object.treatment$family$family!="binomial")){
            stop("Argument \'treatment\' must be a logistic regression (glm object) or a multinomial regression (nnet::multinom)\n",
                 " or a character variable giving the name of the treatment variable. \n")
        }

        if(any(levels(mydata[[treatment]]) %in% mydata[[treatment]] == FALSE)){
            mistreat <- levels(mydata[[treatment]])[levels(mydata[[treatment]]) %in% mydata[[treatment]] == FALSE]
            stop("Argument \'mydata\' needs to take all possible treatment values when using IPTW/AIPTW \n",
                 "missing treatment values: \"",paste(mistreat,collapse="\" \""),"\"\n")
        }

        if(any(is.na(coef(object.treatment)))){
            stop("Cannot handle missing values in the model coefficients (object.treatment) \n")
        }

        if(inherits(object.treatment,"glm") && object.treatment$family$family == "binomial" && length(unique(mydata[[treatment]]))>2){
            stop("Cannot use a logistic regression in argument \"treatment\" when there are more than two treatment categories \n",
                 "Consider subsetting the dataset to have only two treatment categories or using a multinomial regression (nnet::multinom).")
        }
    }

    ## ** object.censor
    if(attr(estimator,"IPCW")){
        ## require censoring model
        if(!inherits(object.censor,"coxph") && !inherits(object.censor,"cph") && !inherits(object.censor,"phreg")){
            stop("Argument \'object.censor\' must be a Cox model \n")
        }
        
        if(inherits(object.event,"glm") && any(n.censor > 0) && all(object.event$prior.weights==1)){
            stop("Argument \'object.event\' should not be a standard logistic regression in presence of censoring \n")
        }

        if(!identical(censorVar.time,eventVar.time)){
            stop("The time variable should be the same in \'object.event\' and \'object.censor\' \n")
        }
        
        if(any(cause == level.censoring)){
            stop("The level indicating censoring should differ between the outcome model and the censoring model \n",
                 "maybe you have forgotten to set the event type == 0 in the censoring model \n")
        }
        
        if(!identical(censorVar.status,eventVar.status)){
            if(sum(diag(table(mydata[[censorVar.status]] == level.censoring, mydata[[eventVar.status]] %in% level.states == FALSE)))!=NROW(mydata)){
                stop("The status variables in object.event and object.censor are inconsistent \n")
            }
        }
        if(any(is.na(coef(object.censor)))){
            stop("Cannot handle missing values in the model coefficients (object.treatment) \n")
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
        ## if(any(data.index != 1:na.omit(n.obs)[1])){
        ##     stop("Argument \'data.index\' cannot be used when using bootstrap resampling \n")
        ## }
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
        if(length(unique(stats::na.omit(n.obs[-1])))>1){
            stop("Arguments \'",paste(names(stats::na.omit(n.obs[-1])), collapse ="\' "),"\' must be fitted using the same number of observations \n")
        }

        if(is.null(data.index)){
            stop("Incompatible number of observations between argument \'data\' and the datasets from argument(s) \'",paste(names(stats::na.omit(n.obs[-1])), collapse ="\' "),"\' \n",
                 "Consider specifying argument \'data.index\' \n \n")
        }
        if(!is.numeric(data.index) || any(duplicated(data.index)) || any(is.na(data.index)) || any(data.index %in% 1:max(n.obs[-1],na.rm = TRUE) == FALSE)){
            stop("Incorrect specification of argument \'data.index\' \n",
                 "Must be a vector of integers between 0 and ",max(n.obs[-1],na.rm = TRUE)," \n")
        }
        
        if(!is.null(object.event)){
            candidateMethods <- paste("predictRiskIID",class(object.event),sep=".")
            if (all(match(candidateMethods,options$method.predictRiskIID,nomatch=0)==0)){
                stop(paste("Could not find predictRiskIID S3-method for ",class(object.event),collapse=" ,"),"\n",
                     "Functional delta method not implemented for this type of object \n",
                     "Set argument \'B\' to a positive integer to use a boostrap instead \n",sep="")
            }
        }
        if(!is.null(object.treatment) && stats::nobs(object.treatment)!=NROW(mydata)){ ## note: only check that the datasets have the same size
            stop("Argument \'treatment\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                 "(discrepancy found in number of rows) \n")
        }

        if(!is.null(object.censor) && coxN(object.censor)!=NROW(mydata)){ ## note: only check that the datasets have the same size
                stop("Argument \'censor\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                     "(discrepancy found in number of rows) \n")
        }

    }

    
    ## ** treatment
    if(!is.null(treatment)){
        if(treatment %in% names(mydata) == FALSE){
            stop("The data set does not seem to have a variable ",treatment," (argument: object.treatment). \n")
        }
        if(is.numeric(mydata[[treatment]])){
            stop("The treatment variable must be a factor variable. \n",
                 "Convert treatment to factor, re-fit the object using this new variable and then call ate. \n")
        }
    }else if(is.null(strata)){
        stop("The treatment variable must be specified using the argument \'treatment\' \n")
    }

    ## ** strata
    if(!is.null(strata)){
        if(attr(estimator, "IPTW")){
            stop("Argument \'strata\' only compatible with the G-formula estimator \n")
        }
        if(any(strata %in% names(mydata) == FALSE)){
            stop("The data set does not seem to have a variable \"",paste0(strata, collapse = "\" \""),"\" (argument: strata). \n")
        }
        if(length(strata) != 1){
            stop("Argument strata should have length 1. \n")
        }
        if(TD == TRUE){
            stop("Landmark analysis is not available when argument strata is specified. \n")
        }
    }
    
    ## ** event time
    if(!is.na(eventVar.time)){
        if(eventVar.time %in% names(mydata) == FALSE){
            stop("The data set does not seem to have a variable ",eventVar.time," (argument: object.event[1]). \n")
        }
        if(is.na(cause)){
            stop("Argument \'cause\' not specified\n")
        }
        data.status <- mydata[[eventVar.status]]
        if(!is.null(treatment)){
            data.strata <- mydata[[treatment]]
        }else{
            data.strata <- mydata[[strata]]
        }
        if(is.factor(data.strata)){
            data.strata <- droplevels(data.strata)
        }
        freq.event <- tapply(data.status, data.strata, function(x){mean(x==cause)})
        count.event <- tapply(data.status, data.strata, function(x){sum(x==cause)})

        if(any(freq.event < 0.01) || any(count.event < 5)  ){
            warning("Rare event, possible violation of the positivity assumption \n")
        }
        if(any(mydata[[eventVar.time]]<0)){
            stop("The time to event variable should only take positive values \n")
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
    
    ## ** iid.nuisance    
    if((return.iid.nuisance == FALSE) & (iid == TRUE) & (attr(estimator, "augmented") == FALSE)){
        warning("Ignoring the uncertainty associated with the estimation of the nuisance parameters may lead to inconsistent standard errors. \n",
                "Consider using a double robust estimator or setting the argument \'known.nuisance\' to TRUE \n")
    }
    

    ## ** output
    return(TRUE)
}


## * helpers
## ** coef.ate
coef.ate <- function(object, ...){
    name.coef <- interaction(object$meanRisk[["treatment"]],object$meanRisk[["time"]])
    value.coef <- object$meanRisk[[paste0("meanRisk.",object$estimator[1])]]
    return(setNames(value.coef, name.coef))
}

## ** vcov.ate
vcov.ate <- function(object, ...){
    if(is.null(object$iid)){
        stop("Missing iid decomposition in object \n",
             "Consider setting the argument \'iid\' to TRUE when calling the ate function \n")
    }
    name.coef <- interaction(object$meanRisk[["treatment"]],object$meanRisk[["time"]])
    iid <- NULL
    for(iT in names(object$iid[[object$estimator[1]]])){ ## iT <- "T0"
        tempo <- object$iid[[object$estimator[1]]][[iT]]
        colnames(tempo) <- interaction(iT,object$time)
        iid <- cbind(iid,tempo)
    }
    return(crossprod(iid[,names(coef(object)),drop=FALSE]))
}

## ** nobs.multinom
nobs.multinom <- function(object,...){
    NROW(object$residuals)
}
##----------------------------------------------------------------------
### ate.R ends here


