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
#' @param max.time [numeric] maximum time of the response of the fitted data.  Only relevant if
#' model \code{response} element has been removed
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
#' of each subject used to fit the object (dim 1),
#' for each subject in newdata (dim 3),
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
#' predCSC.afterT0 <- predict(CSC.fit, newdata = d, cause = 2, times = ttt[ttt>T0], landmark = T0)
#' predCSC.afterT0

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
                                     max.time = NULL,
                                     ...){


    ## ** deal with specific case
    if(type == "survival" && object$surv.type=="survival"){
        predictor.cox <- if(product.limit>0){"predictCoxPL"}else{"predictCox"}
        return(do.call(predictor.cox,
                       args = list(object$models[["OverallSurvival"]], times = times, newdata = newdata, type = "survival",
                                   keep.strata = keep.strata, keep.newdata = keep.newdata,
                                   se = se, band = band, iid = iid, confint = confint, diag = diag,
                                   average.iid = average.iid, store.iid = store.iid)
                       ))
    }
    
    ## ** prepare
    if(missing(newdata)){
        newdata <- eval(object$call$data)
    }else{
        setDT(newdata)
    }

    if (missing(times)) {
        times = object$times
        if (is.null(times)) {
            stop("times must be specified")
        }
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
    if (missing(cause)) {
        cause <- object$theCause
    }
    if (length(cause) > 1){
        stop(paste0("Can only predict one cause. Provided are: ", 
                    paste(cause, collapse = ", "), sep = ""))
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
    type <- match.arg(type, c("absRisk","survival"))
    if(diag && NROW(newdata)!=n.times){
        stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
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
    if (is.null(max.time)) {
        max.time = max(object$response[,"time"])
    }
    valid.times <- times[ times<= max.time] ## prediction times before the event
    if(length(valid.times) == 0){
        if (is.null(eTimes)) {
            stop("eventTimes was removed from model, but no valid times")
        }
        eventTimes <- eTimes[1] ## at least the first event
    }else{
        eventTimes <- eTimes[eTimes <= max(valid.times)] ## jump times before the last prediction time (that is before the last jump)
        if(length(eventTimes) == 0){eventTimes <- eTimes[1]} # at least the first event

    }
    if (is.null(eventTimes)) {
        stop("eventTimes was removed from model - cannot predict")
    }
    eventTimes <- eTimes[eTimes <= max(times)] ## jump times before the last prediction time (that is before the last jump)
    if(length(eventTimes) == 0){eventTimes <- eTimes[1]} # at least the first event
    
    ## order prediction times
    otimes <- order(times)
    ootimes <- order(otimes)
    needOrder <- !identical(1:length(ootimes),ootimes)

    ## ** extract baseline hazard, linear predictor and strata from Cox models
    new.n <- NROW(newdata)
    nEventTimes <- length(eventTimes)
    nCause <- length(causes)
    if(object$surv.type=="survival"){
        nModel <- 2
    }else{
        nModel <- length(causes)
    }
    ls.hazard <- vector(mode = "list", length = nModel)
    ls.cumhazard <- vector(mode = "list", length = nModel)
    M.eXb <- matrix(NA, nrow = new.n, ncol = nModel)
    M.strata.num <- matrix(NA, nrow = new.n, ncol = nModel)
    M.etimes.max <- matrix(NA, nrow = new.n, ncol = nModel)
    ls.infoVar <- setNames(vector(mode = "list", length = nModel), name.model)

    if(length(unlist(coef(object)))==0){
        ## if there is not covariates (only strata) then set the last eventtime to \infty when the last observation is an event
        ls.lastEventTime <- lapply(object$models, function(iM){ ## iM <- object$models[[1]]
            iTempo <- predictCox(iM, times = 0, keep.infoVar = TRUE)
            return(setNames(iTempo$lastEventTime, iTempo$infoVar$strata.levels))
        })
        if(length(unique(lapply(ls.lastEventTime,names)))<=1){
            attr(eventTimes,"etimes.max") <- apply(do.call(rbind,ls.lastEventTime),2,max)
        }
    }

    for(iterC in 1:nModel){ ## iterC <- 1
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
        attr(M.strata.num,paste0("levels",iterC)) <- ls.infoVar[[iterC]]$strata.levels

        ## last event time by strata
        M.etimes.max[,iterC] <- baseline$lastEventTime[M.strata.num[,iterC]+1]
    }

    ## ** compute CIF (aka absolute risk) or event-free survival
    vec.etimes.max <- apply(M.etimes.max,1,max) ## take the max because if not censored for one cause and last event equal to 1 then we have the full curve
    ncores <- parallel::detectCores()
    if(type == "absRisk"){
        outCpp <- predictCIF_cpp(hazard = ls.hazard, 
                                 cumhazard = ls.cumhazard, 
                                 eXb = M.eXb, 
                                 strata = M.strata.num,
                                 newtimes = if(diag){times}else{sort(times)}, 
                                 etimes = eventTimes, 
                                 etimeMax = vec.etimes.max, 
                                 t0 = landmark,
                                 nEventTimes = nEventTimes,
                                 nNewTimes = n.times, 
                                 nData = new.n,
                                 cause = index.cause - 1, 
                                 nCause = nModel,
                                 survtype = (surv.type=="survival"),
                                 productLimit = product.limit>0,
                                 diag = diag,
                                 exportSurv = (se || band || iid || average.iid),
                                 nCores = ncores)

    }else if(type == "survival" && object$surv.type=="hazard"){
        attr(times,"etimes.max") <- attr(eventTimes,"etimes.max")
        
        return(.predictSurv_CSC(object, times = times, newdata = newdata, ls.hazard = ls.hazard, eXb = M.eXb,
                                etimes = eventTimes, etimeMax = vec.etimes.max, strata = M.strata.num,
                                keep.times = keep.times, keep.strata = keep.strata, keep.newdata = keep.newdata,
                                se = se, band = band, iid = iid, 
                                confint = confint, diag = diag, average.iid = average.iid, 
                                store.iid = store.iid, product.limit = product.limit>0))
    }

    ## ** compute standard error for CIF
    if(se || band || iid || average.iid){
        if(!is.na(landmark)){
            stop("standard error for the conditional survival not implemented \n")
        }

        ## design matrix
        new.LPdata <- list()
        for(iCause in 1:nModel){ ## iCause <- 1
            infoVar <- ls.infoVar[[iCause]]
            if(length(infoVar$lpvars) > 0){
                new.LPdata[[iCause]] <- model.matrix(object$models[[iCause]], data = newdata)
            }else{
                new.LPdata[[iCause]] <- matrix(0, ncol = 1, nrow = new.n)
            }  
        }

        ## linear predictors
        nVar.lp <- unlist(lapply(ls.infoVar,function(m){
            length(m$lpvars)
        }))

        Utimes <- sort(unique(times))

        ## Computation of the influence function and/or the standard error
        export <- c("iid"[(iid+band)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE])
        if(!is.null(attr(average.iid,"factor"))){
            if(diag){
                attr(export,"factor") <- attr(average.iid,"factor")
            }else{
                ## re-order columns according to times
                attr(export,"factor") <- lapply(attr(average.iid,"factor"), function(iF){
                    if(NCOL(iF)>1){
                        return(iF[,otimes,drop=FALSE])
                    }else{
                        return(iF)
                    }
                })
            }
        }

        if(product.limit < 0){ ## disregard uncertainty when CIF>1
            check.cif <- outCpp$cif
        }else{  ## usual computation of the uncertainty even when CIF>1
            check.cif <- 0*outCpp$cif
        }

        out.seCSC <- calcSeCSC(object,
                               cif = check.cif,
                               hazard = ls.hazard,
                               cumhazard = ls.cumhazard,
                               survival = outCpp$survival, ## survival at t-
                               object.time = eventTimes,
                               object.maxtime = vec.etimes.max, 
                               eXb = M.eXb,
                               new.LPdata = new.LPdata,
                               new.strata = M.strata.num,                               
                               times = if(diag){times}else{Utimes},
                               ls.infoVar = ls.infoVar,
                               new.n = new.n,
                               cause = index.cause,
                               nCause = nModel,
                               nVar.lp = nVar.lp,
                               surv.type = surv.type,
                               export = export,
                               store.iid = store.iid,
                               diag = diag)

        ootimes2 <- prodlim::sindex(jump.times = Utimes, eval.times = times)
        needOrder2 <- !identical(1:length(ootimes2),ootimes2)
    }
    
    ## ** gather all outputs
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
        c(iI$lpvars.original, iI$stratavars.original)
    })))
    if(keep.newdata[[1]]==TRUE && length(all.covars)>0){
        if (data.table::is.data.table(newdata))
            out$newdata <- newdata[, all.covars, with = FALSE]
        else
            out$newdata <- newdata[, all.covars,drop=FALSE]
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

    ## ** take care of prediction above 1
    if(product.limit<0 && any(out$absRisk>1)){
        index.above1 <- which(out$absRisk>1)
        ## if(iid+band){
        ##     index2.above1 <- which(out$absRisk>1, arr.ind = TRUE)
        ##     for(iObs in 1:NROW(index2.above1)){
        ##         out$absRisk.iid[,index2.above1[iObs,2],index2.above1[iObs,1]] <- 0
        ##     }
        ## }
        ## if(se+band){
        ##     out$absRisk.se[index.above1] <- 0
        ## }        
        out$absRisk[index.above1] <- 1
    }

    ## ** compute confidence interval
    if(confint){
        out <- stats::confint(out)
    }
    if(band[[1]] && se[[1]]==FALSE){
        out["absRisk.se"] <- NULL
    }
    if(band[[1]] && iid[[1]]==FALSE){
        out["absRisk.iid"] <- NULL
    }

    ## ** export
    if(any(na.omit(as.double(out$absRisk))>1) || any(na.omit(as.double(out$absRisk))<0)){
        if(product.limit){
            warning("Estimated risk outside the range [0,1].\n",
                    "Consider setting the argument \'product.limit\' to FALSE. \n")
        }else{
            warning("Estimated risk outside the range [0,1].\n",
                    "Possible cause: incorrect extrapolation, i.e., time and/or covariates used for the prediction differ from those used to fit the Cox models.\n")
        }
    }
    return(out)
}


## * .predictSurv_CSCe
.predictSurv_CSC <- function(object, times, newdata, type, ls.hazard, eXb, etimes, strata, etimeMax,
                             keep.times, keep.strata, keep.newdata,
                             se, band, iid, confint, diag, average.iid, store.iid, product.limit){

    if(!is.logical(diag)){
        stop("Argument \'diag\' must be of type logical. \n")
    }
    if(any(etimes<0)){
        stop("Cannot handle negative event times. \n")
    }
        
    new.n <- NROW(newdata)
    n.times <- length(times)
    n.times2 <- ifelse(diag>0,1,n.times)
    nCause <- length(object$cause)
    n.sample <- NROW(object$response)

    iid.save <- iid
    iid <- (iid || se)
    if(is.null(store.iid)){
        if(is.iidCox(object)){
            store.iid <- object$models[[1]]$iid$store.iid
        }
    }

    ## ** prepare output        
    n.sample <- stats::nobs(object)
    out <- list(survival = NULL, survival.iid = NULL, survival.average.iid = NULL,
                lastEventTime = NA,
                se = se, band = band, type = "survival", diag = diag)
    class(out) <- "predictCox"

    out$survival <- matrix(1, nrow = new.n, ncol = n.times2)

    ## ** prepare hazard
    index.times <- prodlim::sindex(jump.times = etimes, eval.times = times)
    if(!product.limit){
        ls.cumhazard <- lapply(ls.hazard,function(iHazard){
            colCumSum(iHazard)[index.times[index.times>0],,drop=FALSE]
        })
    }

    ## ** get survival
    for(iObs in 1:new.n){ ## iObs <- 1
        if(diag){
            if(times[iObs]>etimeMax[iObs]){
                out$survival[iObs,1] <- NA
            }else if(index.times[iObs]>0){
                if(product.limit){
                    ihazard <- lapply(1:nCause, function(iC){ls.hazard[[iC]][1:index.times[iObs],strata[iObs,iC]+1]*eXb[iObs,iC]})
                    out$survival[iObs,1] <- prod(1-rowSums(do.call(cbind,ihazard)))
                }else{
                    iCumhazard <- sapply(1:nCause, function(iC){ls.cumhazard[[iC]][iObs,strata[iObs,iC]+1]*eXb[iObs,iC]})
                    out$survival[iObs,1] <- exp(-sum(iCumhazard))
                }                
            }
        }else if(diag==0){
            if(product.limit){
                ihazard <- lapply(1:nCause, function(iC){ls.hazard[[iC]][,strata[iObs,iC]+1]*eXb[iObs,iC]})
                out$survival[iObs,index.times>0] <- cumprod(1-rowSums(do.call(cbind,ihazard)))[index.times[index.times>0]]
            }else{
                iCumhazard <- lapply(1:nCause, function(iC){ls.cumhazard[[iC]][,strata[iObs,iC]+1]*eXb[iObs,iC]})
                out$survival[iObs,index.times>0] <- exp(-Reduce("+",iCumhazard))
            }
            if(any(times > etimeMax[iObs])){
                out$survival[iObs,times > etimeMax[iObs]] <- NA
            }
        }
    }

    ## ** update factor with survival
    if(average.iid){
        factor <- attr(average.iid,"factor")
        if(is.null(factor)){
            attr(average.iid,"factor") <- list(-out$survival)
        }else{
            test.list <- !is.list(factor)
            n.factor <- length(factor)
                    
            for(iFactor in 1:n.factor){ ## iFactor <- 1
                ## when only one column and diag = FALSE, use the same weights at all times
                if((diag == FALSE) && (NCOL(factor[[iFactor]])==1) && (NCOL(factor[[iFactor]])==new.n) && (n.times > 1)){
                    factor[[iFactor]] <- matrix(factor[[iFactor]][,1],
                                                nrow = NROW(factor[[iFactor]]),
                                                ncol = n.times, byrow = FALSE)
                }
                ## check dimensions
                if(any(dim(factor[[iFactor]])!=c(new.n, n.times2))){
                    stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ",new.n,",",n.times2," \n")
                }

                factor[[iFactor]] <- -factor[[iFactor]]*out$survival
            }
            attr(average.iid,"factor") <- factor
        }
    }
    
    ## ** compute iid for survival
    if(iid || se || band || average.iid){

        if(iid || se || band){
            iIid <- array(0, dim = c(n.sample, n.times2, new.n))
        }
        if(average.iid){
            
            if(is.null(factor)){
                iAverageIid <- matrix(0, nrow = n.sample, ncol = n.times2)
            }else{
                iAverageIid <- lapply(1:n.factor, function(iF){matrix(0, nrow = n.sample, ncol = n.times2)})
            }
        }

        tsurvival <- t(out$survival)
        
        for(iC in 1:nCause){ ## iC <- 1
            resTempo <- predictCox(object$models[[iC]],
                                   newdata = newdata,
                                   times = times,
                                   diag = diag,
                                   iid = iid || se || band,
                                   se = TRUE,
                                   average.iid = average.iid,
                                   store.iid = store.iid,
                                   type = "cumhazard")

            if(iid || se || band){
                for(iObs in 1:n.sample){
                    iIid[iObs,,] <- iIid[iObs,,] - resTempo$cumhazard.iid[iObs,,] * tsurvival
                }
            }
            if(average.iid){
                if(is.null(factor)){
                    iAverageIid <- iAverageIid + resTempo$cumhazard.average.iid[[1]]
                }else{
                    for(iFactor in 1:n.factor){
                        iAverageIid[[iFactor]] <- iAverageIid[[iFactor]] + resTempo$cumhazard.average.iid[[iFactor]]
                    }
                }
            }
        }

        if(iid || band){
            out$survival.iid <- iIid
        }
        if(se || band){
            out$survival.se <- t(sqrt(apply(iIid^2,2:3,sum)))
        }
        if(average.iid){
            out$survival.average.iid <- iAverageIid
        ## apply(iIid,1:2,mean) - iAverageIid
        }
    }

    ## ** export
    if (keep.times==TRUE){
        out$times <- times
    }
    if (keep.strata[1]==TRUE){
        out$strata <- interaction(unique(lapply(1:nCause,function(iCause){
            factor(strata[,iCause], levels = 0:(length(attr(strata,paste0("levels",iCause)))-1), labels = attr(strata,paste0("levels",iCause)))
        })))
    }
    if( keep.newdata[1]==TRUE){
        out$newdata <- newdata
    }

    ## ** confidence intervals/bands
    if(confint){
        out <- stats::confint(out)
    }
    if(band[1] && se[1]==FALSE){
        out[paste0(type,".se")] <- NULL
    }
    if(band[1] && iid[1]==FALSE){
        out[paste0(type,".iid")] <- NULL
    }
    out$baseline <- FALSE
    return(out)
}
