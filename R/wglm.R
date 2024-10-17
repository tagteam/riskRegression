### wglm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  1 2020 (14:58) 
## Version: 
## Last-Updated: Oct 17 2024 (09:51) 
##           By: Brice Ozenne
##     Update #: 733
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * wglm (documentation)
#' @title Logistic Regression Using IPCW
#' @description Logistic regression over multiple timepoints
#' where right-censoring is handled using inverse probability of censoring weighting (IPCW). 
#'
#' @param formula.event [formula] a formula with a Surv object on the left hand side and the covariates for the logistic regression on the right hand side.
#' @param formula.censor [formula] an optional formula indicating on its right hand side the covariates for the censoring model.
#' @param times [numeric vector] time points at which to model the probability of experiencing an event.
#' @param data [data.frame] dataset containing the time at which the event occured, the type of event, and regressors used to fit the censoring and logistic models.
#' @param cause [character or numeric] the cause of interest. Defaults to the first cause. 
#' @param fitter [character] routine to fit the Cox regression models.
#' @param ties [character] method used to handle ties when using a Cox model (\code{"breslow"} or \code{"efron"}).
#' Ignored if fitter equals to \code{"prodlim"}.
#' @param product.limit [logical] if \code{TRUE} the survival is computed using the product limit estimator.
#' @param store [vector of length 2] Whether prediction should only be computed for unique covariate sets and mapped back to the original dataset (\code{data="minimal"}) and whether the influence function should be stored in a memory efficient way (\code{iid="minimal"}). Otherwise use \code{data="full"} and/or \code{iid="full"}.
#'
#' @details First, a Cox model is fitted (argument formula.censor)
#' and the censoring probabilities are computed relative to each timepoint (argument times) to obtain the censoring weights.
#' Then, for each timepoint, a logistic regression is fitted with the appropriate censoring weights
#' and where the outcome is the indicator of having experience the event of interest (argument cause) at or before the timepoint.
#' 
#' @return an object of class \code{"wglm"}.
#'

## * wglm (examples)
#' @examples
#' library(survival)
#'
#' #### simulate data ####
#' set.seed(10)
#' n <- 250
#' tau <- 1:5
#' d <- sampleData(n, outcome = "competing.risks")
#' dFull <- d[event!=0] ## (artificially) remove censoring
#' dSurv <- d[event!=2] ## (artificially) remove competing risk
#'
#' #### no censoring ####
#' e.wglm <- wglm(Surv(time,event) ~ X1, 
#'                times = tau, data = dFull, product.limit = TRUE)
#' e.wglm ## same as a logistic regression
#'
#' summary(ate(e.wglm, data = dFull, times = tau, treatment = "X1", verbose = FALSE))
#'
#' #### right-censoring ####
#' ## no covariante in the censoring model (independent censoring)
#' eC.wglm <- wglm(Surv(time,event) ~ X1,
#'                times = tau, data = dSurv, product.limit = TRUE)
#' eC.wglm
#'
#' ## with covariates in the censoring model
#' eC2.wglm <- wglm(Surv(time,event) ~ X1 + X8, formula.censor = ~ X1*X8,
#'                  times = tau, data = dSurv)
#' eC2.wglm
#' 
#' #### Competing risks ####
#' ## here Kaplan-Meier as censoring model
#' eCR.wglm <- wglm(Surv(time,event) ~ X1, formula.censor = ~X1,
#'                  times = tau, data = d)
#' eCR.wglm
#' summary(eCR.wglm)
#' eCR.wglm <- wglm(Surv(time,event) ~ X1, formula.censor = ~X1,
#'                  times = tau, data = d)

## * wglm (code)

#' @export
wglm <- function(formula.event, times, data, formula.censor = ~1, cause = NA,
                 fitter = NULL, ties = NULL, product.limit = NULL, store = NULL){
    
    tol <- 1e-12
    mycall <- match.call()

    ## ** check arguments
    ## *** formula.event
    if(inherits(formula.event,"formula")==FALSE){
        stop("Argument \'formula.event' should be a formula. \n")
    }
    varSurv <- SurvResponseVar(formula.event)

    if(!is.null(attr(varSurv$status,"call")) && attr(varSurv$status,"call")!=varSurv$status){
        stop("Mismatch between input and processed Surv element in argument \'formula.event\'. \n",
             "User-input vs. processed: \"",varSurv$status,"\" vs. \"",attr(varSurv$status,"call"),"\".\n",
             "This error may arise when performing operation on the fly in Surv instead of creating a new variable. \n")
    }
    
    ## *** times
    newname <- paste0("XX_",varSurv$status,".",times,"_XX")
    obs.newname <- paste0("XX_observed.",times,"_XX")
    IPCW.newname <- paste0("XX_IPCW.",times,"_XX")
    n.times <- length(times)

    ## *** data
    if(any(newname %in% names(data))){
        stop("Argument \'data\' should not have a column named \"",paste0(newname[newname %in% names(data)],collapse="\" \""),"\"\n",
             "This name is used internally \n.")
    }
    if(any(obs.newname %in% names(data))){
        stop("Argument \'data\' should not have a column named \"",paste0(obs.newname[obs.newname %in% names(data)],collapse="\" \""),"\"\n",
             "This name is used internally \n.")
    }
    if(any(IPCW.newname %in% names(data))){
        stop("Argument \'data\' should not have a column named \"",paste0(IPCW.newname[IPCW.newname %in% names(data)],collapse="\" \""),"\"\n",
             "This name is used internally \n.")
    }
    if((varSurv$time %in% names(data) == FALSE) || (varSurv$status %in% names(data) == FALSE)){
        stop("Mismatch between argument \'formula.event\' and argument \'data\'. \n",
             "Could not find the variable(s) \'",paste(setdiff(c(varSurv$time,varSurv$status), names(data)), collapse = "\', \'"),"\' \n.")
    }
    if(any(is.na(data))){
        warning("Argument \'data\' contains missing values. \n")
    }

    init.Hist <- Hist(time = data[[varSurv$time]], event = data[[varSurv$status]])
    allStates <- c(attr(init.Hist, "cens.code"),attr(init.Hist, "states"))
    code.cens <- attr(init.Hist, "cens.code")
    any.cens <- any((data[[varSurv$status]] == code.cens)*(data[[varSurv$time]]<=max(times)))

    data.time <- data[[varSurv$time]]
    data.status <- data[[varSurv$status]]
    
    ## *** cause
    if(length(cause)>1){
        stop("Argument \'cause\' should have length 1. \n")
    }else if(is.null(cause) || is.na(cause)){
        cause <- attr(init.Hist, "states")[1]
    }else if(cause %in% attr(init.Hist, "states") == FALSE){
        stop("Invalid cause: should be one of \"",paste(attr(init.Hist, "states"), collapse = "\", \""),"\". \n")
    }

    ## *** formula.censor
    if(inherits(formula.censor,"formula")==FALSE){
        stop("Argument \'formula.censor' should be a formula. \n")
    }
    
    if(length(formula.censor)!=2){
        varSurv2 <- SurvResponseVar(formula.censor)
        if(is.null(varSurv2$operator)){
            stop("Argument \'formula.censor\' should not have a left hand side. \n")
        }
        if(varSurv2$time %in% names(data) == FALSE || varSurv2$status %in% names(data) == FALSE){
            stop("Mismatch between argument \'formula.event\' and argument \'data\'. \n",
                 "Could not find the variable(s) \'",paste(setdiff(c(varSurv2$time,varSurv2$status), names(data)), collapse = "\', \'"),"\' \n.")
        }
        if(any(abs(data[[varSurv2$time]]-data.time)>1e-12)){
            stop("Mismatch between the time variable of argument \'formula.event\' and \'formula.censor\'. \n",
                 "Largest discrepancy: ",max(abs(data[[varSurv2$time]]-data.time)),"\n .")
        }
        data.cens <- stats::model.frame(formula.censor, data)$Surv[,2]
        if(all(data.cens + (data[[varSurv$status]] == code.cens) ==1) ){
            stop("Inconsistency between the status variable of argument \'formula.event\' and \'formula.censor\'. \n",
                 "The status variable in argument \'formula.event\' should indicate events while it should indicate censoring in argument \'formula.censor\'. \n",
                 "Typically formula.event=Surv(time,event)~X whereas formula.censor=Surv(time,event==0)~Z. \n")
        }else if(any(data.cens==1 & (data[[varSurv$status]] != code.cens))){
            stop("Inconsistency between the status variable of argument \'formula.event\' and \'formula.censor\'. \n",
                 "Some observations (e.g. ",which(data.cens==1 & (data[[varSurv$status]] != code.cens))[1],") are non-censored for \'formula.event\' while being an event (i.e. censored) in \'formula.censor\'. \n")
        }else if(any(data.cens==0 & (data[[varSurv$status]] == code.cens))){
            stop("Inconsistency between the status variable of argument \'formula.event\' and \'formula.censor\'. \n",
                 "Some observations (e.g. ",which(data.cens==0 & (data[[varSurv$status]] == code.cens))[1],") are censored for \'formula.event\' while not being an event (i.e. non-censored) in \'formula.censor\'. \n")
        }

        ## remove left hand sideo
        formula.censor <- formula(stats::delete.response(terms(formula.censor)))
    }

    ## *** fitter
    if(is.null(fitter)){
        regvar.censor <- all.vars(formula.censor)
        if(length(regvar.censor)>0 && any(grepl("strata",regvar.censor))){
            fitter <- "coxph"
        }else if(length(regvar.censor)>0 && any(grepl("strat",regvar.censor))){
            fitter <- "cph"
        }else if(length(regvar.censor)==0 || all(sapply(regvar.censor, function(iVar){is.factor(data[[iVar]]) || is.character(data[[iVar]]) || is.logical(data[[iVar]])}))){
            fitter <- "prodlim"
        }else{
            fitter <- "coxph"
        }
    }else{
        fitter <- match.arg(fitter,c("coxph","cph","prodlim"))
    }

    if(any.cens == 0){
        fitter <- "prodlim"
        product.limit <- TRUE
    }else if(any.cens && fitter != "prodlim"){
        allTime.cens <- unique(data.time[data.status==code.cens])
        allTime.event <- unique(data.time[data.status==cause])
        M.dist <- as.matrix(stats::dist(c(allTime.cens,allTime.event)))
        diag(M.dist) <- Inf
        Mred.dist <- M.dist[1:length(allTime.cens),length(allTime.cens) + (1:length(allTime.event)),drop=FALSE]
        if(any(Mred.dist==0)){
            if(any(Mred.dist>0 & Mred.dist<5e-5)){
                message("Possible issue due to difference in event times below 5e-5. \n",
                        "Internally adds 1e-6 to censored event times with fitters other than prodlim to break ties. \n")
            }            
        }else{
            ## no censoring happening at the same time as an event
            Mred.dist <- NULL
        }
    }

    ## *** product limit
    if(is.null(product.limit)){
        product.limit <- fitter == "prodlim"
    }

    ## *** ties
    if(is.null(ties)){
        ties <- ifelse(product.limit,"breslow","efron")
    }else{
        ties <- match.arg(ties, c("breslow","efron"))
    }

    ## ** fit censoring model
    ## NOTE: censoring times survival probability should sum up to 1. This is the case when censoring is always posterior to the event, i.e., one should break ties.
    if(any.cens && fitter %in% c("coxph","cph")){
        if(!is.null(Mred.dist)){ ## add extra time to censoring as it occurs after a possible event
            txt.censorSurv <- paste0("Surv(time = ",varSurv$time,"+ 1e-6*(",varSurv$status,"==",code.cens,"), event = ",varSurv$status,"==",code.cens,")",deparse(formula.censor))
        }else{ ## no ambiguity as no ties
            txt.censorSurv <- paste0("Surv(time = ",varSurv$time,", event = ",varSurv$status,"==",code.cens,")",deparse(formula.censor))
        }
        formula.censorSurv <- stats::as.formula(eval(parse(text = txt.censorSurv)))
        if(fitter == "coxph"){
            object.censor <- coxph(formula = formula.censorSurv, data = data, x = TRUE, y = TRUE, ties = ties)
        }else if(fitter == "cph"){
            object.censor <- cph(formula = formula.censorSurv, data = data, x = TRUE, y = TRUE, method = ties)
        }        
    }else if(any.cens && fitter == "prodlim"){
        formula.censorSurv <- stats::as.formula(eval(parse(text = paste0("Hist(time = ",varSurv$time,", event = ",varSurv$status,")",deparse(formula.censor)))))
        object.censor <- prodlim::prodlim(formula.censorSurv, data = data, reverse = TRUE)
        object.censor$method.ties <- ties ## to be read by coxBaseEstimator in predictCox
    }else{
        object.censor <- NULL
    }

    ## ** create glm object
    n.censor <- vector(mode = "numeric", length = n.times)
    out <- list(fit = vector(mode = "list", length = n.times))
    for(iTime in 1:n.times){## iTime <- 1
        iFormula.glm <- update(formula.event, paste0(newname[iTime],"~."))
        
        n.censor[iTime] <- sum((data.time<times[iTime])*(data.status == code.cens))
        data[[newname[iTime]]] <- (data.time<=times[iTime])*(data.status==cause) ## 1 only if event of interest occured before the current timepoint
        data[[obs.newname[iTime]]] <- (1-(data.status == code.cens)*(data.time<times[iTime])) ## 1 if non-censored (i.e. already dead or censored later) and 0 if censored before the current time

        iIndex <- which(data[[obs.newname[iTime]]]>0)
        if(n.censor[iTime]>0){
            iPred <- predictRisk(object.censor, diag = TRUE, newdata = data[iIndex,,drop=FALSE], times = pmin(data.time[iIndex], times[iTime]) - tol,
                                 type = "survival", product.limit = product.limit)[,1]
        }else{
            iPred <- rep(1, length(iIndex))
        }        
        data[[IPCW.newname[[iTime]]]] <- 0
        data[[IPCW.newname[[iTime]]]][iIndex] <- 1/iPred ## 1/P[not-censored] if not censored and 0 otherwise
        
        out$fit[[iTime]] <- suppressWarnings(do.call(stats::glm, list(formula = iFormula.glm, family = stats::binomial(link = "logit"),
                                                                      data = data, weights = data[[IPCW.newname[[iTime]]]])))
        ## Warning message:
        ##             In eval(family$initialize) : non-integer #successes in a binomial glm!
        out$fit[[iTime]]$time.prior.weights <- pmin(data.time, times[iTime]) - tol
        out$fit[[iTime]]$prior.weights2 <- rep(0, NROW(data))
        out$fit[[iTime]]$prior.weights2[iIndex] <- 1/iPred^2
    }

    ## ** export
    out$call <- mycall
    out$formula <- formula.event
    out$n <- NROW(data)
    out$n.censor <- n.censor
    out$model.censor <- object.censor
    out$data <- data
    out$var.outcome <- varSurv$status
    out$var.time <- varSurv$time
    out$name.IPCW <- IPCW.newname
    out$name.outcome <- newname
    out$times <- times
    out$theCause <- cause
    out$causes <- allStates
    out$product.limit <- product.limit
    out$fitter <- fitter
    class(out) <- append("wglm",class(out))
    return(out)
}

## * nobs.wglm
#' @export
nobs.wglm <- function(object, ...){
    return(object$n)
}
## * formula.wglm
#' @export
formula.wglm <- function(x, ...){
    return(x$formula)
}

## * coef.wglm
#' @title Estimates from IPCW Logistic Regressions
#' @description Display the estimated regression parameters from logistic regressions.
#'
#' @param object a wglm object.
#' @param times [numeric vector] time points at which the estimates should be output. 
#' @param simplify [logical] should the ouput be converted to a vector when only one timepoint is requested. Otherwise will always return a matrix.
#' @param ... Not used.
#' @export
coef.wglm <- function(object, times = NULL, simplify = TRUE, ...){
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% object$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    M.coef <- do.call(rbind,setNames(lapply(object$fit, coef),object$time))
    if(simplify && length(times)==1){
        out <- stats::setNames(M.coef[object$times %in% times,], colnames(M.coef))
    }else{
        out <- M.coef[object$times %in% times,,drop=FALSE]
    }
    return(out)
}

## * vcov.wglm
#' @title Variance-covariance for IPCW Logistic Regressions
#' @description Compute the variance-covariance matrix of the estimated model parameters of IPCW logistic regressions.
#'
#' @param object a wglm object.
#' @param times [numeric vector] time points at which the variance-covariance matrix should be output. 
#' @param simplify [logical] should the ouput be converted to a matrix when only one timepoint is requested. Otherwise will always return a list.
#' @param ... Not used.
#' @export
vcov.wglm <- function(object, times = NULL, simplify = TRUE, ...){
    ls.iid <- iid(object, times = times, simplify = simplify)
    if(is.list(ls.iid)){
        out <- lapply(ls.iid, crossprod)
    }else{
        out <- crossprod(ls.iid)
    }
    return(out)
}

## * summary.wglm
#' @export
summary.wglm <- function(object, print = TRUE, se = "robust", times = NULL, ...){

    se <- match.arg(se, c("robust","model-wknown","robust-wknown"))
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% object$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    if(se == "robust"){
        object.iid <- lava::iid(object, simplify = FALSE, times = times)
    }else if(se == "robust-wknown"){
        object.iid <- lapply(object$fit[match(times, object$times)],lava::iid)
    }
    
    out <- setNames(vector(mode = "list", length = n.times), times)
    if(print){
        print(object, short = TRUE)
    }
    for(iTime in 1:n.times){
        iTime2 <- which(object$times == times[iTime])
        suppressWarnings(out[[iTime]] <- summary(object$fit[[iTime2]])$coef)
        if(se %in% c("robust","robust-wknown")){
            out[[iTime]][,"Std. Error"] <- sqrt(diag(crossprod(object.iid[[iTime]])))
            out[[iTime]][,3] <- out[[iTime]][,"Estimate"]/out[[iTime]][,"Std. Error"] ## name change from z to t stat when using quasibinomial instead binomial
            out[[iTime]][,4] <- 2*(1-pnorm(abs(out[[iTime]][,3])))
        }
        if(print){
            cat("  ------ time: ",times[iTime]," ----------------------------------------------\n",sep="")
            iIPCW <- object$data[[object$name.IPCW[iTime]]]
            cat("  - number (events, no event, censoring): ",sum(object$fit[[iTime]]$y),", ",object$n-sum(object$fit[[iTime]]$y)-object$n.censor[iTime],", ",object$n.censor[iTime],"\n",sep="")
            cat("  - IPCW (min,median,max): ",paste(round(quantile(iIPCW[iIPCW>0], c(0,0.5,1)), digits = 5), collapse = ", "),"\n",sep="")
            cat("\n")
            iPrint <- out[[iTime]]
            rownames(iPrint) <- paste0("  ",rownames(iPrint))
            print(iPrint)
            cat("\n")
        }
    }

    return(invisible(out))
}

## * print.wglm
#' @export
print.wglm <- function(x, times = NULL, short = FALSE, ...){

    ## ** prepare
    if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    
    if(length(x$causes)>1){
        text.CR <- paste0(" for cause ",x$theCause)
    }else{
        text.CR <- ""
    }
    if(any(x$n.censor>0)){
        text.IPCW <- "IPCW "
    }else{
        text.IPCW <- ""
    }

    if(any(x$n.censor>1)){
        ff.censor <- stats::delete.response(stats::terms(formula(x$model.censor)))
    }
    ff.outcome <- stats::delete.response(stats::terms(formula(x$formula)))
    if(!short){
        x.coef <- coef(x, simplify = FALSE)
        M.print <- cbind(n.censor = x$n.censor,
                         n.event = sapply(x$fit, function(iM){sum(iM$y)}),
                         "IPCW(max)" = sapply(x$name.IPCW, function(iName){max(x$data[[iName]])}),
                         x.coef)

        colnames(M.print)[-(1:3)] <- colnames(x.coef)
        rownames(M.print) <- paste0("   ",x$times)    
    }

    ## ** print
    cat("\t\t",text.IPCW,"logistic regression",text.CR," \n\n",sep="")
    cat(" - structure: ",deparse(x$formula[[2]])," with possible states: ",paste(x$causes,collapse=", "),". \n",sep="")
    if(any(x$n.censor>1)){
        cat(" - censoring model: ",deparse(ff.censor)," (fitter: ",x$fitter," function). \n",sep="")
    }
    cat(" - outcome model: ",deparse(ff.outcome)," (fitter: glm function). \n",sep="")
    cat(" - estimated regression parameters: \n")
    if(!short){
        print(M.print)
    }
    cat("\n")

    return(NULL)
}

## * score.wglm
#' @title Score for IPCW Logistic Regressions
#' @description Compute the first derivative of the log-likelihood for IPCW logistic regressions.
#'
#' @param x a wglm object.
#' @param indiv [logical] should the individual score be output? Otherwise the total score (i.e. summed over all individuals will be output).
#' @param times [numeric vector] time points at which the score should be output. 
#' @param simplify [logical] should the ouput be converted to a matrix when only one timepoint is requested. Otherwise will always return a list.
#' @param ... Not used.
#' @export
score.wglm <- function(x, indiv = FALSE, times = NULL, simplify = TRUE, ...){
    if(inherits(x,"glm") && !inherits(x,"wglm")){
        x <- list(fit = list("1" = x),
                  times = "1")
        times <- "1"
    }else if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        X <- stats::model.matrix(iObject)
        pi <- stats::predict(iObject, type = "response")
        Y <- iObject$y
        W <- iObject$prior.weights
        out[[iTime]] <- colMultiply_cpp(X, scale = W*(Y - pi))
        colnames(out[[iTime]]) <- colnames(X)
        if(indiv==FALSE){out[[iTime]] <- colSums(out[[iTime]])}
    }

    if(n.times == 1 && simplify){
        return(out[[1]])
    }else{
        return(out)
    }
}

## * information.wglm
#' @title Information for IPCW Logistic Regressions
#' @description Compute the information (i.e. opposite of the expectation of the second derivative of the log-likelihood) for IPCW logistic regressions.
#'
#' @param x a wglm object.
#' @param times [numeric vector] time points at which the score should be output. 
#' @param simplify [logical] should the ouput be converted to a matrix when only one timepoint is requested. Otherwise will always return a list.
#' @param ... Not used.
#' @export
information.wglm <- function(x, times = NULL, simplify = TRUE, ...){
    if(inherits(x,"glm") && !inherits(x,"wglm")){
        x <- list(fit = list("1" = x),
                  times = "1")
        times <- "1"
    }else if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        X <- stats::model.matrix(iObject)
        pi <- stats::predict(iObject, type = "response")
        W <- iObject$prior.weights
        out[[iTime]] <- t(colMultiply_cpp(X, scale = W*pi*(1-pi))) %*% X
        rownames(out[[iTime]]) <- colnames(out[[iTime]])
    }

    if(n.times == 1 && simplify){
        return(out[[1]])
    }else{
        return(out)
    }

}

## * iid.wglm
#' @title IID for IPCW Logistic Regressions
#' @description Compute the decomposition in iid elements of the ML estimor of IPCW logistic regressions.
#'
#' @param x a wglm object.
#' @param times [numeric vector] time points at which the iid should be output. 
#' @param simplify [logical] should the ouput be converted to a matrix when only one timepoint is requested. Otherwise will always return a list.
#' @param ... Not used.
#' @export
iid.wglm <- function(x, times = NULL, simplify = TRUE, ...){
    if(inherits(x,"glm") && !inherits(x,"wglm")){
        x <- list(fit = list("1" = x),
                  times = "1",
                  n.censor = 0)
        times <- "1"
        class(x) <- append("wglm",class(x))
    }else if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    ls.score <- lava::score(x, times = times, simplify = FALSE, indiv = TRUE)
    ls.info <- lava::information(x, times = times, simplify = FALSE)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        ## ** compute the uncertainty related to the weights
        ## S score, I information matrix, H hessian, X design matrix, Y outcome, \pi fitted probabilities
        ## \beta regresssion parameters (logit) \eta nuisance parameters (censoring)
        ## The estimating equation of \beta is
        ## S = 0 = sum_i W_i X_i (Y_i - \pi_i(\beta))
        ##       = n\int W(O,\Prob_n) X(O) (Y(O) - \pi(O,\beta(\Prob_n))) d\Prob_n(O)
        ## dS(\Prob_h)/dh =   n\int dW(O, \Prob_h)/dh X(O) (Y(O) - \pi(O,\beta(\Prob_h))) d\Prob_h(O)
        ##                  - n\int W(O, \eta(\Prob_h)) X(O) d\pi(O,\beta(\Prob_h))/dh d\Prob_h(O)
        ##                  + n\int W(O, \eta(\Prob_h)) X(O) (Y(O) - \pi(O,\beta(\Prob_h))) d(d\Prob_h(O)/dh)
        ## dS(\Prob_h)/dh =   n(\int dW(O,\Prob_h)/dh X(O) (Y(O) - \pi(O)) d\Prob_h(O))
        ##                  - n(\int W(0)) X(O) d\pi(0,\beta)/d\beta d\Prob_h(O)) * (d\beta(\Prob_h)/d\Prob_h)
        ##                  + n\int W(0) X(O) (Y(O) - \pi(O)) d(d\Prob_n(O)/dh)
        ## \Prob_h = (1-h)\Prob + h \delta_{O_i}
        ## dS(\Prob_h)/dh|h = 0 = n\int IF_{W_O}(O_i) X(O) (Y(O) - \pi(O)) d\Prob_h(O) 
        ##                        + dS/d\beta IF_\beta(O_i)
        ##                        + (-S + S_i) where S=0 (since it is at ML)
        ## IF_\beta(O_i) = (S_i + n\int IF_{W_O}(O_i) X(O) (Y(O) - \pi(O)) d\Prob_h(O)) / (-dS/d\beta)
        ##               = (S_i + n*AIF_{W,X*(Y-pi)}(O_i) / I
        ## NOTE: IF_{W_O}(O_i) = dW(O,\eta)/d\eta \IF_\eta(O_i)
        
        if(x$n.censor[iTime]>0){
            n.obs <- sum(coxN(iObject))
            W2 <- iObject$prior.weights2
            X <- stats::model.matrix(iObject)
            pi <- stats::predict(iObject, type = "response")
            Y <- iObject$y
            factor <- TRUE        

            attr(factor,"factor") <- lapply(apply(colMultiply_cpp(X, -W2*(Y - pi)), 2, list), function(iVec){cbind(iVec[[1]])})
            iPred <- predictRisk(x$model.censor, diag = TRUE, newdata = x$data, times = iObject$time.prior.weights,
                                 type = "survival", product.limit = x$product.limit, average.iid = factor)
        }

        ## ** assemble uncertainty
        ## (S+dS/dW)/I
        if(x$n.censor[iTime]>0){
            out[[iTime]] <- (ls.score[[iTime]] + do.call(cbind,attr(iPred, "average.iid"))*NROW(x$data)) %*% solve(ls.info[[iTime]])
        }else{
            out[[iTime]] <- ls.score[[iTime]] %*% solve(ls.info[[iTime]])
        }
        ## ## WRONG VERSION
        ## factor1 <- colMultiply_cpp(X, (Y - pi)) %*% iVcov
        ## factor2 <- - ls.score[[iTime]] %*% crossprod(iVcov) %*% t(colMultiply_cpp(X, scale = -pi*(1-pi))) %*% X
        
        ## factor <- TRUE        
        ## attr(factor,"factor") <- lapply(apply(factor1 + factor2, 2, list), function(iVec){cbind(iVec[[1]])})
        ## iPred <- predictRisk(x$model.censor, diag = TRUE, newdata = x$data, times = iObject$time.prior.weights,
        ##                      type = "survival", product.limit = FALSE, average.iid = factor)
        ## out[[iTime]] <- out[[iTime]] - do.call(cbind,attr(iPred, "average.iid"))*NROW(x$data)

    }
    
    ## ** export
    if(n.times == 1 && simplify){
        return(out[[1]])
    }else{
        return(out)
    }
}


## * weights.wglm
#' @title Extract IPCW Weights
#' @description Extract IPCW weights of IPCW logistic regressions.
#'
#' @param object a wglm object.
#' @param times [numeric vector] time points at which the weights should be output. 
#' @param simplify [logical] should the ouput be converted to a vector when only one timepoint is requested. Otherwise will always return a matrix.
#' @param ... Not used.
#' @export
#'
weights.wglm <- function(object, times = NULL, simplify = TRUE, ...){

    out <- object$data[object$name.IPCW]
    if(!is.null(times)){
        if(any(times %in% object$time == FALSE)){
            stop("Unknown timepoint ",paste(setdiff(times,object$time), collapse = ", ")," in argument \'times\'. \n",
                 "Valid timepoints: ",paste(object$time, collapse = ", "),". \n")
        }
        out <- out[match(times,object$time)]
    }

    if(simplify && NCOL(out)==1){
        return(out[[1]])
    }else{
        return(out)
    }
}
######################################################################
### wglm.R ends here
