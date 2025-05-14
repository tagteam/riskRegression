### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: May 14 2025 (16:10) 
##           By: Thomas Alexander Gerds
##     Update #: 2580
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
#' @description Use the g-formula or the IPW or the double robust estimator
#' to estimate the average treatment
#'     effect (absolute risk difference or ratio)
#' based on Cox regression with or without competing risks.
#' @name ate
#' 
#' @param event Outcome model which describes how the probability of experiencing a terminal event depends
#'     on treatment and covariates. The object carry its own call and have a \code{predictRisk} method. See details section and examples section.
#' @param treatment Name of the treatment variable or treatment model which describes how the probability of being allocated to a treatment group depends
#'     on covariates. See details section and examples section.
#' @param censor Censoring model which describes how the probability of being censored depends
#'     on treatment and covariates. See details section and examples section.
#' @param data [data.frame or data.table] Data set in which to evaluate risk predictions
#' based on the outcome model
#' @param data.index [numeric vector] Position of the observation in argument data relative to the dataset used to obtain the argument event, treatment, censor.
#' Only necessary for the standard errors when computing the Average Treatment Effects on a subset of the data set. 
#' @param contrasts [character vector] levels of the treatment variable for which the risks should be assessed and compared. Default is to consider all levels.
#' @param allContrasts [2-row character matrix] levels of the treatment variable to be compared. Default is to consider all pairwise comparisons.
#' @param strata [character] Strata variable on which to compute the average risk.
#' Incompatible with treatment. Experimental.
#' 
#' @param times [numeric vector] Time points at which to evaluate average treatment effects.
#' @param cause [integer/character] the cause of interest.
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
#' \code{"foreach"} is the default and the only option on Windows. It uses \code{parallel} to create a cluster.
#' Other operating systems can use \code{"mclapply"}.
#' This argument is ignored when \code{mc.cores=1} and \code{cl=NULL}.
#' @param mc.cores [integer, >0] The number of cores to use,
#' i.e., the upper limit for the number of child processes that run simultaneously.
#' Passed to \code{parallel::mclapply} or \code{parallel::makeCluster}.
#' The option is initialized from environment variable mc_cores if set.
#' @param cl A parallel socket cluster used to perform cluster calculation in parallel (output by \code{parallel::makeCluster}).
#' The packages necessary to run the computations (e.g. riskRegression) must already be loaded on each worker.
#' Only used when \code{handler="foreach"}.
#' @param verbose [logical] If \code{TRUE} inform about estimated run time.
#' @param ... passed to predictRisk
#'
#' @author Brice Ozenne \email{broz@@sund.ku.dk}
#' and Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#'
#' @details Argument \bold{event}: \itemize{
#' \item IPTW estimator: a character vector indicating the event and status variables.
#' \item G-formula or AIPTW estimator: a survival model (e.g. Cox \code{"survival::coxph"}, Kaplan Meier \code{"prodlim::prodlim"}),
#' a competing risk model (e.g. cause-specific Cox \code{"riskRegression::CSC"}),
#' a logistic model (e.g. \code{"stats::glm"} when argument \code{censor} is \code{NULL} otherwise \code{"riskRegression::wglm"}).
#' Other models can be specified, provided a suitable \code{predictRisk} method exists, but the standard error should be computed using non-parametric bootstrap,
#' as the influence function of the estimator will typically not be available.
#' }
#' Argument \bold{treatment}: \itemize{
#' \item G-formula estimator: a character indicating the treatment variable.
#' \item IPTW or AIPTW estimator: a \code{"stats::glm"} model with family \code{"binomial"} (two treatment options) or a \code{"nnet::multinom"} (more than two treatment options).
#' }
#' Argument \bold{censor}: \itemize{
#' \item G-formula estimator: NULL
#' \item IPTW or AIPTW estimator: NULL if no censoring and otherwise a survival model (e.g. Cox \code{"survival::coxph"}, Kaplan Meier \code{"prodlim::prodlim"}) 
#' }
#' 
#' Argument \bold{estimator}: when set to \code{"AIPTW"} with argument \code{event} being IPCW logistic model (\code{"riskRegression::wglm"}), the integral term w.r.t. to the martingale of the censoring process is not computed,
#' i.e. using the notation of Ozenne et al. (2020) an AIPTW,IPCW estimator instead of an AIPTW,AIPCW is evaluated. \cr
#'
#' In presence of censoring, the computation time and memory usage for the evaluation of the AIPTW estimator and its uncertainty do not scale well with the number of observations (n) or the number of unique timepoints (T).
#' Point estimation involves n by T matrices, influence function involves n by T by n arrays. \itemize{
#' \item for large datasets (e.g. n>5000), bootstrap is recommended as the memory need for the influence function is likely prohibitive.
#' \item it is possible to decrease the memory usage for the point estimation by setting the (hidden) argument \code{store=c(size.split=1000)}. The integral term of the AIPTW estimator is then evaluated for 1000 observations at a time, i.e. involing matrices of size 1000 by T instead of n by T. This may lead to increased computation time.
#' \item reducing the number of unique timepoints (e.g. by rounding them) will lead to an approximate estimation procedure that is less demanding both in term of computation and memory. The resulting estimator will be more variable than the one based on the original timepoints (i.e. wider confidence intervals).
#' }
#'  
#' 
#' 
#' 
#' 
#' @references
#' Brice Maxime Hugues Ozenne, Thomas Harder Scheike, Laila Staerk, and Thomas
#'    Alexander Gerds. On the estimation of average treatment effects with right-
#'    censored time to event outcome and competing risks. Biometrical Journal, 62
#'    (3):751--763, 2020.
#' 
#' @seealso
#' \code{\link{autoplot.ate}} for a graphical representation the standardized risks. \cr
#' \code{\link{coef.ate}} to output estimates for the the average risk, or difference in average risks, or ratio between average risks. \cr
#' \code{\link{confint.ate}} to output a list containing all estimates (average & difference & ratio) with their confidence intervals and p-values. \cr
#' \code{\link{model.tables.ate}} to output a data.frame containing one type of estimates (average or difference or ratio) with its confidence intervals and p-values.   \cr
#' \code{\link{summary.ate}} for displaying in the console a summary of the results.

## * ate (examples)
#' @examples 
#' library(survival)
#' library(prodlim)
#' library(data.table)
#'
#' ############################
#' #### Survival settings  ####
#' #### ATE with Cox model ####
#' ############################
#'
#' #### generate data ####
#' n <- 100
#' set.seed(10)
#' dtS <- sampleData(n, outcome="survival")
#' dtS$time <- round(dtS$time,1)
#' dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
#'
#' ##### estimate the survival model ####
#' ## Cox proportional hazard model
#' fit.Cox <- coxph(Surv(time,event)~ X1+X2, data=dtS, y=TRUE, x=TRUE)
#' ## Kaplan Meier - same as copxh(Surv(time,event)~ strata(X1)+strata(X2), ties = "breslow")
#' fit.KM <- prodlim(Hist(time,event)~ X1+X2, data=dtS)
#'
#' ##### compute the ATE ####
#' ## at times 5, 6, 7, and 8 using X1 as the treatment variable
#' ## standard error computed using the influence function
#' ## confidence intervals / p-values based on asymptotic results
#' ateFit1a <- ate(fit.Cox, treatment = "X1", times = 5:8, data = dtS)
#' summary(ateFit1a)
#' model.tables(ateFit1a, type = "meanRisk")
#' model.tables(ateFit1a, type = "diffRisk")
#' model.tables(ateFit1a, type = "ratioRisk")
#'
#' ## relaxing PH assumption using a stratified model
#' ateFit1a.NPH <- ate(fit.KM, treatment = "X1", times = 5:8, data = dtS,
#'                     product.limit = FALSE)
#' model.tables(ateFit1a.NPH, times = 5) ## no PH assumption
#' model.tables(ateFit1a, times = 5)  ## PH assumption
#' 
#' \dontrun{
#' #### ATE with confidence bands ####
#' ## same as before with in addition the confidence bands / adjusted p-values
#' ## (argument band = TRUE)
#' ateFit1b <- ate(fit.Cox, treatment = "X1", times = 5:8, data = dtS, band = TRUE)
#' summary(ateFit1b)
#'
#' ## by default bands/adjuste p-values computed separately for each treatment modality
#' summary(ateFit1b, band = 1,
#'         se = FALSE, type = "diffRisk", short = TRUE, quantile = TRUE)
#' ## adjustment over treatment and time using the band argument of confint
#' summary(ateFit1b, band = 2,
#'        se = FALSE, type = "diffRisk", short = TRUE, quantile = TRUE)
#' }
#' 
#' \dontrun{
#' #### ATE with non-parametric bootstrap ####
#' ## confidence intervals / p-values computed using 1000 bootstrap samples
#' ## (argument se = TRUE and B = 1000) 
#' ateFit1c <- ate(fit.Cox, treatment = "X1", times = 5:8, data = dtS,
#'                 se = TRUE, B = 50, handler = "mclapply")
#' ## NOTE: for real applications 50 bootstrap samples is not enough 
#'
#' ## same but using 2 cpus for generating and analyzing the bootstrap samples
#' ## (parallel computation, argument mc.cores = 2) 
#' ateFit1d <- ate(fit.Cox, treatment = "X1", times = 5:8, data = dtS, 
#'                 se = TRUE, B = 50, mc.cores = 2)
#'
#' ## manually defining the cluster to be used
#' ## useful when specific packages need to be loaded in each cluster
#' fit.CoxNL <- coxph(Surv(time,event)~ X1+X2+rcs(X6),data=dtS,y=TRUE,x=TRUE)
#'
#' cl <- parallel::makeCluster(2)
#' parallel::clusterEvalQ(cl, library(rms))
#' 
#' ateFit1e <- ate(fit.CoxNL, treatment = "X1", times = 5:8, data = dtS, 
#'                 se = TRUE, B = 50, handler = "foreach", cl = cl)
#' parallel::stopCluster(cl)
#' }
#'
#'
#' ################################################
#' #### Competing risks settings               ####
#' #### ATE with cause specific Cox regression ####
#' ################################################
#' 
#' #### generate data ####
#' n <- 500
#' set.seed(10)
#' dt <- sampleData(n, outcome="competing.risks")
#' dt$X1 <- factor(rbinom(n, prob = c(0.2,0.3) , size = 2), labels = paste0("T",0:2))
#'
#' #### estimate cause specific Cox model ####
#' fit.CR <-  CSC(Hist(time,event)~ X1+X8,data=dt,cause=1)
#' 
#' #### compute the ATE ####
#' ## at times 1, 5, 10
#' ## using X1 as the treatment variable
#' ateFit2a <- ate(fit.CR, treatment = "X1", times = c(1,5,10), data = dt, 
#'                 cause = 1, se = TRUE, band = TRUE)
#' summary(ateFit2a)
#' data.table::as.data.table(ateFit2a)
#' 
#' #############################################
#' #### Survival settings without censoring ####
#' #### ATE with glm                        ####
#' #############################################
#' 
#' #### generate data ####
#' n <- 100
#' dtB <- sampleData(n, outcome="binary")
#' dtB[, X2 := as.numeric(X2)]
#' 
#' ##### estimate a logistic regression model ####
#' fit.glm <- glm(formula = Y ~ X1+X2, data=dtB, family = "binomial")
#' 
#' #### compute the ATE ####
#' ## using X1 as the treatment variable
#' ## only point estimate (argument se = FALSE)
#' ateFit1a <- ate(fit.glm, treatment = "X1", data = dtB, se = FALSE)
#' ateFit1a
#' 
#' \dontrun{
#' ## with confidence intervals
#' ateFit1b <- ate(fit.glm, data = dtB, treatment = "X1",
#'                times = 5) ## just for having a nice output not used in computations
#' summary(ateFit1b, short = TRUE)
#' 
#' ## using the lava package
#' library(lava)
#' ateLava <- estimate(fit.glm, function(p, data){
#' a <- p["(Intercept)"] ; b <- p["X11"] ; c <- p["X2"] ;
#' R.X11 <- expit(a + b + c * data[["X2"]])
#' R.X10 <- expit(a + c * data[["X2"]])
#' list(risk0=R.X10,risk1=R.X11,riskdiff=R.X11-R.X10)},
#' average=TRUE)
#' ateLava
#' }
#'
#' ############################
#' #### Survival settings  ####
#' #### ATE with glm       ####
#' ############################
#' 
#' ## see wglm for handling right-censoring with glm
#'
#' #################################
#' #### Double robust estimator ####
#' #################################
#' 
#' \dontrun{
#' ## generate data
#' n <- 500
#' set.seed(10)
#' dt <- sampleData(n, outcome="competing.risks")
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
#'                  cause = 1)
#' summary(ateRobust)
#' 
#' ## compare various estimators
#' system.time( ## about 1.5s
#' ateRobust2 <- ate(event = m.event,
#'                  treatment = m.treatment,
#'                  censor = m.censor,
#'                  estimator = c("GFORMULA","IPTW","AIPTW"),
#'                  data = dt, times = c(5:10), 
#'                  cause = 1, se = TRUE)
#' )
#' data.table::as.data.table(ateRobust2, type = "meanRisk")
#' data.table::as.data.table(ateRobust2, type = "diffRisk")
#'
#' ## reduce memory load by computing the integral term
#' ## on only 100 observations at a time (only relevant when iid = FALSE)
#' ateRobust3 <- ate(event = m.event,
#'                  treatment = m.treatment,
#'                  censor = m.censor,
#'                  data = dt, times = 5:10, 
#'                  cause = 1, se = FALSE,
#'                  store = c("size.split" = 100), verbose = 2)
#' coef(ateRobust3) - coef(ateRobust) ## same
#'
#' ## approximation to speed up calculations
#' dt$time.round <- round(dt$time/0.5)*0.5 ## round to the nearest half
#' dt$time.round[dt$time.round==0] <- 0.1 ## ensure strictly positive event times
#' mRound.event <-  CSC(Hist(time.round,event)~ X1+X2+X3+X5+X8,data=dt)
#' mRound.censor <-  coxph(Surv(time.round,event==0)~ X1+X2+X3+X5+X8,data=dt, x = TRUE, y = TRUE)
#' system.time( ## about 0.4s
#' ateRobustFast <- ate(event = mRound.event, treatment = m.treatment,
#'                     censor = mRound.censor,
#'                  estimator = c("GFORMULA","IPTW","AIPTW"), product.limit = FALSE,
#'                  data = dt, times = c(5:10), cause = 1, se = TRUE)
#' )
#' coef(ateRobustFast) - coef(ateRobust) ## inaccuracy
#' }
#' 
## * ate (code)
#' @rdname ate
#' @export
ate <- function(event,
                treatment,
                censor = NULL,
                data,
                data.index = NULL,
                estimator = NULL,
                strata = NULL,
                contrasts = NULL,
                allContrasts = NULL,
                times,
                cause = NA,
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
    call <- match.call()

    ## ** initialize arguments
    if(missing(data)){ 
        if(inherits(event,"formula") || inherits(treatment,"treatment") || inherits(censor,"censor")){ ## Support to mice: try to find variables in the parent environment
            etc.vars <- c(all.vars(event),
                          all.vars(treatment),
                          all.vars(censor))
            if(all(etc.vars %in% ls(envir = parent.frame()))){
                etc.formula <- paste("~",paste(etc.vars, collapse ="+"))
                etc.mf <- parse(text=paste0("stats::model.frame(",etc.formula,")"))
                data <- data.table::as.data.table(eval(etc.mf, envir = parent.frame()))
            }else{
                stop("Argument \"data\" is missing and could not find variables \"",paste(setdiff(etc.vars,ls(envir = parent.frame())), collapse = "\", \""),"\" in parent frame.\n ",
                     "Consider specifying the argument \"data\". \n")
            }
        }else{
            stop("Argument \"data\" is missing, with no default. \n")
        }
    }else if(is.data.table(data)){
        data <- data.table::copy(data)
    }else{
        data <- data.table::as.data.table(data)
    }

    init <- ate_initArgs(object.event = event,
                         object.treatment = treatment,
                         object.censor = censor,
                         mydata = data,
                         data.index = data.index,
                         estimator = estimator,
                         times = times,
                         cause = cause,
                         handler = handler,
                         product.limit = dots$product.limit,
                         store = dots$store)

    object.event <- init$object.event
    object.treatment <- init$object.treatment
    object.censor <- init$object.censor
    data.index <- init$data.index
    times <- init$times
    handler <- init$handler
    treatment <- init$treatment
    cause <- init$cause
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
    store <- init$store

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
    check <- ate_checkArgs(call = call,
                           object.event = object.event,
                           object.censor = object.censor,
                           object.treatment = object.treatment,
                           mydata = data,
                           treatment = treatment,
                           strata = strata,
                           contrasts = contrasts,
                           allContrasts = allContrasts,
                           times = times,
                           cause = cause,
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
    if(is.null(allContrasts)){
        allContrasts <- utils::combn(contrasts, m = 2)
    }else if(any(contrasts %in% allContrasts == FALSE)){
        contrasts <- contrasts[contrasts %in% allContrasts]
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
                allContrasts = allContrasts,
                computation.time = list("point" = NA, "iid" = NA, "bootstrap" = NA), 
                inference = data.frame("se" = se, "iid" = iid, "band" = band, "bootstrap" = B>0, "ci" = FALSE, "p.value" = FALSE,
                                       "conf.level" = NA, "alternative" = NA,
                                       "bootci.method" = NA, "n.bootstrap" = if(B>0){B}else{NA},
                                       "method.band" = NA, "n.sim" = NA)
                )
    attr(out$theCause,"cause") <- level.states
    attr(out$theCause,"level.censoring") <- level.censoring
    if(!is.na(eventVar.time)){
        attr(out$variables,"range.time") <- range(data[[eventVar.time]])

        if(!is.null(strata)){
            var.group <- strata
        }else{
            var.group <- treatment
        }
        data.contrasts <- data[data[[var.group]] %in% contrasts]

        ## do not use [,,by=...] because there can be confusion between the argument times
        # and a column named times when present in the data.table
        data.timeContrasts <- do.call(rbind,by(data.contrasts, INDICES = data.contrasts[[var.group]], FUN = function(iDF){
            iM <- do.call(rbind,lapply(times, function(t){c(times = t, pc = sum(iDF[[eventVar.time]]>=t))}))
            return(cbind(iDF[1,.SD, .SDcols = var.group], iM))                                                               
        }))
    attr(out$eval.times,"n.at.risk") <- dcast(data.timeContrasts,
                                              value.var = "pc",
                                              formula = as.formula(paste0(var.group,"~times")))
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
            object.event <- iidCox(object.event, tau.max = max(times), store.iid = store[["iid"]], return.object = TRUE)
        }

        if(attr(estimator,"IPCW")){
            if(identical(is.iidCox(object.censor),FALSE)){
                if(verbose>1/2){cat(" censoring")}
                object.censor <- iidCox(object.censor, store.iid = store[["iid"]], tau.max = max(times))
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
                                 allContrasts=allContrasts,
                                 levels=levels,
                                 times=times,
                                 cause=cause,
                                 n.censor = n.censor,
                                 level.censoring = level.censoring,
                                 estimator = estimator,
                                 eventVar.time = eventVar.time,
                                 eventVar.status = eventVar.status,
                                 censorVar.time = censorVar.time,
                                 censorVar.status = censorVar.status,
                                 return.iid.nuisance = return.iid.nuisance,
                                 data.index = data.index,
                                 method.iid = method.iid,
                                 store = store,
                                 verbose = verbose),
                            dots[names(dots) %in% "store" == FALSE])
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
        
        key1 <- c("treatment","time")
        key2 <- c("treatment.A","treatment.B","time")

        if(B>0){
                                        # {{{ Bootstrap
            ## *** Bootstrap
            ## display start
            if (verbose>1/2){ 
                cat(" - Non-parametric bootstrap using ",B," samples and ",mc.cores," core",if(mc.cores>1){"s"},"\n", sep ="")
                cat("                            (expected time: ",round(as.numeric(out$computation.time$point)*B/mc.cores,2)," seconds)\n", sep = "")
            }

            name.estimate <- c(paste("mean",out$meanRisk$estimator,out$meanRisk$time,out$meanRisk$treatment,sep="."),
                               paste("difference",out$diffRisk$estimator,out$diffRisk$time,out$diffRisk$A,out$diffRisk$B,sep="."),
                               paste("ratio",out$ratioRisk$estimator,out$ratioRisk$time,out$ratioRisk$A,out$ratioRisk$B,sep="."))
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
                                   cl = cl)

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
                                IPTW = pointEstimate$store$iid.IPTW,
                                AIPTW = pointEstimate$store$iid.AIPTW)
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
        out[c("meanRisk","diffRisk","ratioRisk","inference","inference.allContrasts","inference.contrasts","transform")] <- stats::confint(out, p.value = TRUE)

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
##----------------------------------------------------------------------
### ate.R ends here


