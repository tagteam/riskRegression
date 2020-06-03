## * predict.CauseSpecificCox (documentation)
#' @title Predicting Absolute Risk from Cause-Specific Cox Models
#' @description  Apply formula to combine two or more Cox models into absolute risk (cumulative incidence function).
#' @name predict.CauseSpecificCox
#' @aliases predict.CauseSpecificCox
#' @aliases predictBig.CauseSpecificCox
#' 
#' @param object The fitted cause specific Cox model
#' @param newdata [data.frame or data.table]  Contain the values of the predictor variables
#' defining subject specific predictions relative to each cause.
#' Should have the same structure as the data set used to fit the \code{object}.
#' @param times [numeric vector] Time points at which to return
#' the estimated absolute risk.
#' @param cause [integer/character] Identifies the cause of interest among the competing
#'     events.
#' @param type [character] Can be changed to \code{"survival"} if the event free survival should be output instead of the absolute risk.
#' @param landmark [integer] The starting time for the computation of the cumulative risk.
#' @param keep.times [logical] If \code{TRUE} add the evaluation times to the output.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the covariates used to make the prediction in the output list. 
#' @param keep.strata [logical] If \code{TRUE} add the value of the strata used to make the prediction in the output list. 
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param average.iid [logical]. If \code{TRUE} add the average of the influence function over \code{newdata} to the output.
#' @param product.limit [logical]. If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param store.iid [character] Implementation used to estimate the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}. 
#' @param diag [logical] when \code{FALSE} the absolute risk/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param ... not used.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds
#'     tag@@biostat.ku.dk
#' 
#' @details
#' This function computes the absolute risk as given by formula 2 of (Ozenne et al., 2017).
#' Confidence intervals and confidence bands can be computed using a first order von Mises expansion.
#' See the section "Construction of the confidence intervals" in (Ozenne et al., 2017).
#' 
#' A detailed explanation about the meaning of the argument \code{store.iid} can be found
#' in (Ozenne et al., 2017) Appendix B "Saving the influence functions".
#' 
#' Note: for Cox regression models with time varying
#'     covariates it does not make sense to use this function, because
#'     the predicted risk has to be a measurable function of the data
#'     available at the time origin.
#' 
#' The iid decomposition is output using an array containing the value of the influence
#' of each subject used to fit the object (dim 3),
#' for each subject in newdata (dim 1),
#' and each time (dim 2).
#' 
#' @seealso
#' \code{\link{confint.predictCSC}} to compute confidence intervals/bands.
#' \code{\link{autoplot.predictCSC}} to display the predictions.
#' 
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#'

## * predict.CauseSpecificCox (examples)
#' @examples
#' library(survival)
##' library(prodlim)
#' #### generate data ####
#' set.seed(5)
#' d <- sampleData(80,outcome="comp") ## training dataset
#' nd <- sampleData(4,outcome="comp") ## validation dataset
#' d$time <- round(d$time,1) ## create tied events
#' ttt <- sort(sample(x = unique(d$time), size = 10))
#'
#' ## estimate a CSC model based on the coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X3+X8, data=d, method = "breslow")
#'
#' ## compute the absolute risk of cause 1, in the validation dataset
#' ## at time 1:10
#' CSC.risk <-  predict(CSC.fit, newdata=nd, times=1:10, cause=1)
#' CSC.risk
#'
#' ## compute absolute risks with CI for cause 2
#' ## (without displaying the value of the covariates)
#' predict(CSC.fit,newdata=nd,times=1:10,cause=2,se=TRUE,
#'         keep.newdata = FALSE)
#'
#' ## other example
#' library(survival)
#' CSC.fit.s <- CSC(list(Hist(time,event)~ strata(X1)+X2+X9,
#'  Hist(time,event)~ X2+strata(X4)+X8+X7),data=d, method = "breslow")
#' predict(CSC.fit.s,cause=1,times=ttt,se=1L) ## note: absRisk>1 due to small number of observations
#' 
#' ## using the cph function instead of coxph
#' CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", fitter = "cph")#' 
#' predict(CSC.cph, newdata = d, cause = 2, times = ttt)
#' 
#' ## landmark analysis
#' T0 <- 1
#' predCSC_afterT0 <- predict(CSC.fit, newdata = d, cause = 2, times = ttt[ttt>T0], landmark = T0)
#' predCSC_afterT0

## * predict.CauseSpecificCox (code)
#' @rdname predict.CauseSpecificCox
#' @method predict CauseSpecificCox
#' @export
predict.CauseSpecificCox <- function(object,
                                     newdata,
                                     times,
                                     cause,
                                     type = "absRisk",
                                     landmark = NA,
                                     keep.times = 1L,
                                     keep.newdata = 1L,
                                     keep.strata = 1L,
                                     se = FALSE,
                                     band = FALSE,
                                     iid = FALSE,
                                     confint = (se+band)>0,
                                     average.iid = FALSE,
                                     product.limit = TRUE,
                                     store.iid = "full",
                                     diag = FALSE,
                                     ...){


    if(missing(newdata)){
        newdata <- eval(object$call$data)
    }else{
        setDT(newdata)
    }

    ## ** event-free survival instead of absolute risk
    type <- match.arg(type, c("absRisk","survival"))
    if(type=="survival"){
        return(.predictSurv_CSC(object, newdata = newdata, times = times, product.limit = product.limit, se = se, diag = diag, iid = iid, average.iid = average.iid))
    }

    ## ** prepare
    n.times <- length(times)
    if(object$fitter=="phreg"){newdata$entry <- 0} 
    new.n <- NROW(newdata)
    ## if(data.table::is.data.table(newdata)){
    ## newdata <- data.table::copy(newdata)
    ## }else{
    ## newdata <- data.table::as.data.table(newdata)
    ## }
    
    surv.type <- object$surv.type
    if (length(cause) > 1){
        stop(paste0("Can only predict one cause. Provided are: ", 
                    paste(cause, collapse = ", "), sep = ""))
    }
    if (missing(cause)) {
        cause <- object$theCause
    }
	
    ## causes
    # NOTE: cannot use only eventtimes of cause 1 otherwise wrong estimation of the survival in the absolute risk
    causes <- object$causes
    index.cause <- which(causes == cause)
    name.model <- names(object$models)
    
    ## event times
    eTimes <- object$eventTimes

    ## ** check
    if (any(match(as.character(cause), causes, nomatch = 0)==0L))
        stop(paste0("Cannot find all requested cause(s) ...\n\n", 
                    "Requested cause(s): ", paste0(cause, collapse = ", "), 
                    "\n Available causes: ", paste(causes, collapse = ", "), 
                    "\n"))
    if (surv.type == "survival") {
        if (object$theCause != cause) 
            stop("Object can be used to predict cause ", object$theCause, 
                 " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
    }
    if(any(is.na(times))){
        stop("NA values in argument \'times\' \n")
    }
    if(length(landmark)!=1){
        stop("\'t0\' must have length one \n")
    }
    if(!is.logical(diag)){ 
        stop("Argument \'diag\' must be logical \n")
    }
    
    if(diag){
        if(iid[[1]]==TRUE && store.iid[[1]] == "minimal"){
            stop("Arguments \'store.iid\' must equal \"full\" when \'diag\' is TRUE \n")
        }
        if(NROW(newdata)!=n.times){
            stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
        }
    }
    if(average.iid==TRUE && !is.null(attr(average.iid,"factor"))){
        if(iid && !is.null(attr(average.iid,"factor"))){
            stop("Attribute \"factor\" of argument \'average.iid\' not available when \'iid\' is TRUE \n")
        }
        if(se && !is.null(attr(average.iid,"factor"))){
            stop("Attribute \"factor\" of argument \'average.iid\' not available when \'se\' is TRUE \n")
        }

        test.list <- !is.list(attr(average.iid,"factor"))
        if(test.list){
            stop("Attribute \"factor\" of argument \'average.iid\' must be a list \n")
        }
        test.matrix <- any(unlist(lapply(attr(average.iid,"factor"), is.matrix))==FALSE)
        if(test.matrix){
            stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices \n")
        }
        for(iFactor in 1:length(attr(average.iid,"factor"))){ ## iFactor <- 1
            ## check dimensions
            if(NROW(attr(average.iid,"factor")[[iFactor]])!=new.n){
                stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices with ",new.n," rows \n")
            }
            if(NCOL(attr(average.iid,"factor")[[iFactor]]) %in% c(1, diag + (1-diag)*n.times) == FALSE){
                stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices with ",diag + (1-diag)*n.times," columns\n")
            }
        }
    }
    
    ## relevant event times to use
    valid.times <- times[times<=max(object$response[,"time"])] ## prediction times before the event 
    if(length(valid.times) == 0){
        eventTimes <- eTimes[1] ## at least the first event
    }else{
        eventTimes <- eTimes[eTimes <= max(valid.times)] ## jump times before the last prediction time (that is before the last jump)
        if(length(eventTimes) == 0){eventTimes <- eTimes[1]} # at least the first event
    }
    
    ## order prediction times
    otimes <- order(times)
    ootimes <- order(otimes)
    needOrder <- !identical(1:length(ootimes),ootimes)

    ## ** extract baseline hazard, linear predictor and strata from Cox models
    new.n <- NROW(newdata)
    nEventTimes <- length(eventTimes)
    nCause <- length(causes)
        
    ls.hazard <- vector(mode = "list", length = nCause)
    ls.cumhazard <- vector(mode = "list", length = nCause)
    M.eXb <- matrix(NA, nrow = new.n, ncol = nCause)
    M.strata.num <- matrix(NA, nrow = new.n, ncol = nCause)
    M.etimes.max <- matrix(NA, nrow = new.n, ncol = nCause)
    ls.infoVar <- setNames(vector(mode = "list", length = nCause), name.model)
     
    for(iterC in 1:nCause){ ## iterC <- 1

        ## when surv.type = "hazard" and iterC corresponds to the cause and no se/iid
        ## we could only compute cumhazard (i.e. not compute hazard).
        ## But since computing hazard has little impact on the performance it is done anyway
        baseline <- predictCox(object$models[[iterC]],
                               centered = FALSE,
                               times = eventTimes,
                               newdata = NULL,
                               type = c("hazard","cumhazard"), 
                               keep.strata = TRUE,
                               keep.times = TRUE,
                               se = FALSE,
                               keep.infoVar = TRUE)
        ls.infoVar[[iterC]] <- baseline$infoVar

        ## baseline hazard from the Cox model
        ls.cumhazard[[iterC]] <- matrix(baseline$cumhazard, byrow = FALSE, nrow = nEventTimes)
        ls.hazard[[iterC]] <- matrix(baseline$hazard, byrow = FALSE, nrow = nEventTimes)
          
        ## linear predictor for the new observations
        M.eXb[,iterC] <- exp(coxLP(object$models[[iterC]], data = newdata, center = FALSE))

        ## strata for the new observations
        strataTempo <- coxStrata(object$models[[iterC]], data = newdata, 
                                 sterms = ls.infoVar[[iterC]]$strata.sterms, 
                                 strata.vars = ls.infoVar[[iterC]]$stratavars, 
                                 strata.levels = ls.infoVar[[iterC]]$strata.levels)
        M.strata.num[,iterC] <- as.numeric(strataTempo) - 1
        ## last time by strata
        M.etimes.max[,iterC] <- baseline$lastEventTime[M.strata.num[,iterC]+1]
    }

    ## ** compute CIF
    vec.etimes.max <- apply(M.etimes.max,1,min)
    outCpp <- predictCIF_cpp(hazard = ls.hazard, 
                             cumhazard = ls.cumhazard, 
                             eXb = M.eXb, 
                             strata = M.strata.num,
                             newtimes = if(diag){times}else{sort(times)}, 
                             etimes = eventTimes, 
                             etimeMax = vec.etimes.max, 
                             t0 = landmark,
                             nEventTimes = nEventTimes,
                             nNewTimes = length(times), 
                             nData = new.n,
                             cause = index.cause - 1, 
                             nCause = nCause,
                             survtype = (surv.type=="survival"),
                             productLimit = product.limit,
                             diag = diag,
                             exportSurv = (se || band || iid || average.iid))

    ## ** compute standard error for CIF
    if(se || band || iid || average.iid){
        if(!is.na(landmark)){
            stop("standard error for the conditional survival not implemented \n")
        }

        ## design matrix
        new.LPdata <- list()
        for(iCause in 1:nCause){ ## iCause <- 1
            infoVar <- ls.infoVar[[iCause]]
            if(length(infoVar$lpvars) > 0){
                new.LPdata[[iCause]] <- model.matrix(object$models[[iCause]], data = newdata)
            }else{
                new.LPdata[[iCause]] <- matrix(0, ncol = 1, nrow = new.n)
            }  
        }

        ## linear predictors
        nVar <- unlist(lapply(ls.infoVar,function(m){
            length(m$lpvars)
        }))

        Utimes <- sort(unique(times))

        ## influence function for the parameters of the Cox model
        if(is.iidCox(object)){
            for(iModel in 1:nCause){ # iModel <- 1
                if(!identical(object$models[[iModel]]$iid$time,eventTimes)){
                    object$models[[iModel]]$iid <- selectJump(object$models[[iModel]]$iid, times = eventTimes,
                                                              type = c("hazard","cumhazard"))
                    store.iid <- unique(c(store.iid,object$models[[iModel]]$iid$store.iid))
                }
            }
            if(length(store.iid)!=1){
                stop("When pre-computed, the iid for the Cox models must be done using the same \'store.iid\' argument \n")
            }
        }else{
            object <- iidCox(object, tau.hazard = eventTimes,
                             store.iid = store.iid, return.object = TRUE)
        }

        ## Computations
        if(se || band || iid){
            out.seCSC <- calcSeCSC(object,
                                   cif = outCpp$cif,
                                   hazard = ls.hazard,
                                   cumhazard = ls.cumhazard,
                                   survival = outCpp$survival,
                                   object.time = eventTimes,
                                   object.maxtime = vec.etimes.max, 
                                   eXb = M.eXb,
                                   new.LPdata = new.LPdata,
                                   new.strata = M.strata.num,                               
                                   times = if(diag){times}else{Utimes},
                                   ls.infoVar = ls.infoVar,
                                   new.n = new.n,
                                   cause = index.cause,
                                   nCause = nCause,
                                   nVar = nVar,
                                   surv.type = surv.type,
                                   export = c("iid"[(iid+band)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE]),
                                   store.iid = store.iid,
                                   diag = diag)
        }else if(average.iid){ ## only average.iid in export
            if(!is.null(attr(average.iid,"factor"))){
                if(diag){
                    factor <- attr(average.iid,"factor")
                }else{
                    ## re-order columns according to times
                    factor <- lapply(attr(average.iid,"factor"), function(iF){
                        if(NCOL(iF)>1){return(iF[,otimes,drop=FALSE])}else{return(iF)}
                    })
                }
            }
            out.seCSC <- calcAiidCSC(object,
                                     cif = outCpp$cif,
                                     hazard = ls.hazard,
                                     cumhazard = ls.cumhazard,
                                     survival = outCpp$survival,
                                     object.time = eventTimes,
                                     object.maxtime = vec.etimes.max, 
                                     eXb = M.eXb,
                                     new.LPdata = new.LPdata,
                                     new.strata = M.strata.num,                               
                                     times = if(diag){times}else{Utimes},
                                     ls.infoVar = ls.infoVar,
                                     new.n = new.n,
                                     cause = index.cause,
                                     nCause = nCause,
                                     nVar = nVar,
                                     surv.type = surv.type,
                                     factor = factor,
                                     store.iid = store.iid,
                                     diag = diag)
        }

        ootimes2 <- prodlim::sindex(jump.times = Utimes, eval.times = times)
        needOrder2 <- !identical(1:length(ootimes2),ootimes2)
    }
    
    ## ** export
    ## gather all outputs
    if(needOrder && (diag == FALSE)){
        out <- list(absRisk = outCpp$cif[,ootimes,drop=FALSE]) # reorder prediction times
    }else{
        out <- list(absRisk = outCpp$cif) # reorder prediction times
    }
    
    if(se+band){
        if(needOrder2 && (diag == FALSE)){
            out$absRisk.se <- out.seCSC$se[,ootimes2,drop=FALSE]
        }else{
            out$absRisk.se <- out.seCSC$se
        }
    }
    if(iid+band){
        if(needOrder2 && (diag == FALSE)){
            out$absRisk.iid <- out.seCSC$iid[,ootimes2,,drop=FALSE]
        }else{
            out$absRisk.iid <- out.seCSC$iid
        }
    }
    if(average.iid){
        if(needOrder2 && (diag == FALSE)){
            if(is.list(out.seCSC$average.iid)){
                out$absRisk.average.iid <- lapply(out.seCSC$average.iid, function(iIID){iIID[,ootimes2,drop=FALSE]})
            }else{
                out$absRisk.average.iid <- out.seCSC$average.iid[,ootimes2,drop=FALSE]
            }
        }else{
            out$absRisk.average.iid <- out.seCSC$average.iid
        }
        if(is.list(attr(average.iid,"factor"))){
            names(out$absRisk.average.iid) <- names(attr(average.iid,"factor"))
        }
    }
    if(keep.times){out$times <- times}

    all.covars <- unique(unlist(lapply(ls.infoVar, function(iI){
        c(iI$lpvars.original, iI$strata.vars.original)
    })))
    if(keep.newdata[[1]]==TRUE && length(all.covars)>0){
        out$newdata <- newdata[, all.covars, with = FALSE]
    }
    if(keep.strata==TRUE){
        allStrata <- unique(unlist(lapply(ls.infoVar,"[[","stratavars.original")))
        if (length(allStrata)>0){
            newdata <- copy(newdata[,allStrata, with = FALSE])
            newdata[, (allStrata) := lapply(allStrata, function(col){paste0(col,"=",.SD[[col]])})]
            out$strata <- newdata[, interaction(.SD, sep = " "), .SDcols = allStrata]
        }
    }
    out$diag <- diag
    out$se <- se
    out$keep.times <- keep.times
    out$band <- band
    
    class(out) <- "predictCSC"

    ## compute confidence interval
    if(confint){
        out <- stats::confint(out)
    }
    if(band[[1]] && se[[1]]==FALSE){
        out["absRisk.se"] <- NULL
    }
    if(band[[1]] && iid[[1]]==FALSE){
        out["absRisk.iid"] <- NULL
    }
    ## export
    return(out)
}





## * .predictSurv_CSC
.predictSurv_CSC <- function(object, newdata, times, product.limit, diag, iid, average.iid, se){

    if(object$surv.type=="survival"){
        ## names(object$models)
        predictor.cox <- if(product.limit){"predictCoxPL"}else{"predictCox"}

        out <- do.call(predictor.cox,
                       args = list(object$models[["OverallSurvival"]], newdata = newdata, times = times, diag = diag, iid = iid, average.iid = average.iid, type = "survival", se = se)
                       )
        
    }else if(object$surv.type=="hazard"){
        if(!is.logical(diag)){
            stop("Argument \'diag\' must be of type logical \n")
        }
        
        n.obs <- NROW(newdata)
        n.times <- length(times)
        n.cause <- length(object$cause)
        n.data <- NROW(object$response)

        iid.save <- iid
        iid <- (iid || se)

        ## ** prepare output        
        out <- list(survival = NULL, survival.iid = NULL, survival.average.iid = NULL,
                    lastEventTime = NA,
                    se = FALSE, band = FALSE, type = "survival", diag = diag, times = times)
        class(out) <- "predictCox"

        if(diag){
            out$survival <- matrix(0, nrow = n.obs, ncol = 1)
            if(se){out$survival.se <- matrix(0, nrow = n.obs, ncol = 1)}
            if(iid){out$survival.iid <- array(0, dim = c(n.obs, 1, n.data))}
            if(average.iid){out$survival.average.iid <- matrix(0, nrow = n.data, ncol = 1)}
        }else{
            out$survival <- matrix(0, nrow = n.obs, ncol = n.times)
            if(se){out$survival.se <- matrix(0, nrow = n.obs, ncol = n.times)}
            if(iid){out$survival.iid <- array(0, dim = c(n.obs, n.times, n.data))}
            if(average.iid){out$survival.average.iid <- matrix(0, nrow = n.data, ncol = n.times)}
        }

        ## ** compute survival && iid for cumulative hazard
        if(product.limit){

            ## find all jumps
            jump.time <- object$eventTime[object$eventTime <= max(times)]
            if(0 %in% jump.time){
                jumpA.time <- c(jump.time,max(object$eventTime)+1e-10)
            }else{
                jumpA.time <- c(0,jump.time,max(object$eventTime)+1e-10)
            }
            n.jumpA <- length(jumpA.time)

            ## compute hazard at all times
            predAll.hazard <- matrix(0, nrow = n.obs, ncol = n.jumpA)
            for(iC in 1:n.cause){
                outHazard <- predictCox(object$models[[iC]],
                                        newdata = newdata,
                                        times = jump.time,
                                        iid = FALSE, se = FALSE,
                                        type = "hazard")

                if(0 %in% jump.time){
                    predAll.hazard <- predAll.hazard + cbind(outHazard$hazard,NA)
                }else{
                    predAll.hazard <- predAll.hazard + cbind(0,outHazard$hazard,NA)
                }

                if(iid){
                    outHazard2 <- predictCox(object$models[[iC]],
                                             newdata = newdata,
                                             times = times,
                                             iid = iid,
                                             diag = diag,
                                             type = "cumhazard")
                    out$survival.iid <- out$survival.iid + outHazard2$cumhazard.iid
                }
            }

            index.jump <- prodlim::sindex(eval.times = times,
                                          jump.times = jumpA.time)

            out$survival <- rowCumProd(1-predAll.hazard)[,index.jump,drop=FALSE]
            if(diag){
                out$survival <- cbind(diag(out$survival))
            }
            
        }else{

            for(iC in 1:n.cause){
                outHazard <- predictCox(object$models[[iC]],
                                        newdata = newdata,
                                        times = times,
                                        diag = diag,
                                        iid = iid,
                                        type = "cumhazard")
                out$survival <- out$survival + outHazard$cumhazard
                if(iid){
                    out$survival.iid <- out$survival.iid + outHazard$cumhazard.iid
                }
            }
            out$survival <- exp(-out$survival)
        }

        ## ** compute iid for survival
        if(iid){
            if(diag){
                out$survival.iid[,1,] <- colMultiply_cpp(-out$survival.iid[,1,], scale = out$survival)
            }else{
                for(iTau in 1:n.times){ ## iTau <- 1
                    if(n.obs==1){
                        out$survival.iid[,iTau,] <- -out$survival.iid[,iTau,] * out$survival[,iTau]
                    }else{
                        out$survival.iid[,iTau,] <- colMultiply_cpp(-out$survival.iid[,iTau,], scale = out$survival[,iTau])
                    }
                }            
            }
        }

        ## ** compute se for survival
        if(se){
            out$survival.se[] <- sqrt(apply(out$survival.iid^2,1:2,sum))
            if(iid.save == FALSE){
                out$survival.iid <- NULL
            }
        }

        ## ** compute average iid survival
        if(average.iid){
            if(iid && !is.null(attr(average.iid,"factor"))){
                stop("Attribute \"factor\" of argument \'average.iid\' not available when \'iid\' is TRUE \n")
            }
            if(se && !is.null(attr(average.iid,"factor"))){
                stop("Attribute \"factor\" of argument \'average.iid\' not available when \'se\' is TRUE \n")
            }

                    
            if(iid){
                out$survival.average.iid <- t(apply(out$survival.iid,2:3,mean))
            }else{
                factor <- attr(average.iid,"factor")
                if(is.null(factor)){
                    factor <- list(-out$survival)
                    n.factor <- 1
                }else{
                    test.list <- !is.list(factor)
                    if(test.list){
                        stop("Attribute \"factor\" of argument \'average.iid\' must be a list \n")
                    }
                    test.matrix <- any(unlist(lapply(factor, is.matrix))==FALSE)
                    if(test.matrix){
                        stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices \n")
                    }
                    n.factor <- length(factor)
                    out$survival.average.iid <- vector(mode = "list", length = n.factor)
                    
                    for(iFactor in 1:n.factor){ ## iFactor <- 1
                        ## when only one column and diag = FALSE, use the same weights at all times
                        if((diag == FALSE) && (NCOL(factor[[iFactor]])==1) && (NCOL(factor[[iFactor]])==n.obs) && (n.times > 1)){
                            factor[[iFactor]] <- matrix(factor[[iFactor]][,1],
                                                        nrow = NROW(factor[[iFactor]]),
                                                        ncol = n.times, byrow = FALSE)
                        }
                        ## check dimensions
                        if(any(dim(factor[[iFactor]])!=c(n.obs, diag + (1-diag)*n.times))){
                            stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ",n.obs,",",diag + (1-diag)*n.times," \n")
                        }

                        factor[[iFactor]] <- -factor[[iFactor]]*out$survival
                        out$survival.average.iid[[iFactor]] <- matrix(0, nrow = n.data, ncol = n.times)
                    }
                    ## stop("factor for average.iid when calling predict.CauseSpecificCox with type=\"survival\" not implemented \n")
                }
                attr(average.iid,"factor") <- factor
                
                for(iC in 1:n.cause){
                    iAvIID <- predictCox(object$models[[iC]],
                                         newdata = newdata,
                                         times = times,
                                         diag = diag,
                                         iid = FALSE,
                                         average.iid = average.iid,
                                         type = "cumhazard")$cumhazard.average.iid

                    if(is.list(out$survival.average.iid)){
                        for(iF in 1:n.factor){ ## iF <- 1
                            out$survival.average.iid[[iF]] <- out$survival.average.iid[[iF]] + iAvIID[[iF]]
                        }
                    }else{
                        out$survival.average.iid <- out$survival.average.iid + iAvIID[[1]]
                    }
                }
            }
        }
        
    }

    ## ** export
    return(out)
}
