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
#' @param landmark [integer] The starting time for the computation of the cumulative risk.,
#' @param keep.times [logical] If \code{TRUE} add the evaluation times to the output.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the covariates used to make the prediction in the output list. 
#' @param keep.strata [logical] If \code{TRUE} add the value of the strata used to make the prediction in the output list. 
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param average.iid [logical]. If \code{TRUE} add the average of the influence function over \code{newdata} to the output.
#' @param product.limit [logical]. If true the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param store.iid [character] Implementation used to estimate the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}. 
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
#' @rdname predict.CauseSpecificCox
#' @examples
#' library(survival)
#' 
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
                                     ...){
    if(object$fitter=="phreg"){newdata$entry <- 0} 
    if(missing(newdata)){newdata <- eval(object$call$data)}
    data.table::setDT(newdata)
    
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
    
    ## event times
    eTimes <- object$eventTimes
    
    if (any(match(as.character(cause), causes, nomatch = 0)==0L))
        stop(paste0("Cannot find all requested cause(s) ...\n\n", 
                    "Requested cause(s): ", paste0(cause, collapse = ", "), 
                    "\n Available causes: ", paste(causes, collapse = ", "), 
                    "\n"))
    ## stopifnot(match(as.character(cause), causes, nomatch = 0) != 
    ## 0)
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
    if(any(is.na(unlist(coef(object))))){
        stop("Incorrect object",
             "One or several model parameters have been estimated to be NA \n")
    }
  
    # relevant event times to use  
    eventTimes <- eTimes[which(eTimes <= max(times))] 
    if(length(eventTimes) == 0){eventTimes <- min(times)} # at least the first event

    # order prediction times
    ootimes <- order(order(times))

                                        # predict cumulative cause specific hazards
    new.n <- NROW(newdata)
    nEventTimes <- length(eventTimes)
    nCause <- length(causes)
        
    if (surv.type == "hazard") {
        
        ls.hazard <- vector(mode = "list", length = nCause)
        ls.cumhazard <- vector(mode = "list", length = nCause)
        M.eXb <- matrix(NA, nrow = new.n, ncol = nCause)
        M.strata <- matrix(NA, nrow = new.n, ncol = nCause)
        M.etimes.max <- matrix(NA, nrow = new.n, ncol = nCause)
        ls.infoVar <- setNames(vector(mode = "list", length = nCause),as.character(causes))
     
        for(iterC in 1:nCause){

            if(iterC == index.cause || product.limit || se || iid || average.iid){            
                typeC <- c("hazard","cumhazard")
            }else{
                typeC <- "cumhazard"
            }
            baseline <- predictCox(object$models[[iterC]], centered = FALSE,
                                   times = eventTimes, newdata = NULL,
                                   type = typeC, 
                                   keep.strata = TRUE, keep.times = TRUE,
                                   se = FALSE, keep.infoVar = TRUE)
            ls.infoVar[[iterC]] <- baseline$infoVar

            ## baseline hazard from the Cox model
            ls.cumhazard[[iterC]] <- matrix(baseline$cumhazard, byrow = FALSE, nrow = nEventTimes)
            if("hazard" %in% typeC){            
                ls.hazard[[iterC]] <- matrix(baseline$hazard, byrow = FALSE, nrow = nEventTimes)
            }else{
                ls.hazard[[iterC]] <- matrix()
            }
          
            ## linear predictor for the new observations
            M.eXb[,iterC] <- exp(coxLP(object$models[[iterC]], data = newdata, center = FALSE))
          
            ## strata for the new observations
            M.strata[,iterC] <- as.numeric(coxStrata(object$models[[iterC]], data = newdata, 
                                                     sterms = ls.infoVar[[iterC]]$sterms, 
                                                     strata.vars = ls.infoVar[[iterC]]$strata.vars, 
                                                     levels = levels(baseline$strata), 
                                                     strata.levels = ls.infoVar[[iterC]]$strata.levels))-1
          
            ## last time by strata
            M.etimes.max[,iterC] <- baseline$lastEventTime[M.strata[,iterC]+1]
        }
        
        
    }else if (surv.type == "survival"){
        ls.infoVar <- setNames(list(NULL, NULL), c(as.character(cause), "OverallSurvival"))
        #### cause ####
        
        ## baseline hazard from the Cox model
        baseline_Cause <- predictCox(object$models[[paste("Cause",cause)]],
                                     centered = FALSE,
                                     times = eventTimes,
                                     newdata = NULL,
                                     type = c("hazard","cumhazard"),
                                     keep.strata = TRUE,
                                     keep.times = TRUE,
                                     se = FALSE, keep.infoVar = TRUE)
        ls.infoVar[[as.character(cause)]] <- baseline_Cause$infoVar
        
        ## linear predictor for the new observations
        eXb_Cause <- cbind(exp(coxLP(object$models[[paste("Cause",cause)]], data = newdata, center = FALSE)))
        
        ## strata for the new observations
        strata_Cause <- coxStrata(object$models[[paste("Cause",cause)]],
                                  data = newdata,
                                  sterms = ls.infoVar[[as.character(cause)]]$infoVar$sterms,
                                  strata.vars = ls.infoVar[[as.character(cause)]]$infoVar$strata.vars,
                                  levels = levels(baseline_Cause$strata),
                                  strata.levels = ls.infoVar[[as.character(cause)]]$infoVar$strata.levels)
        
        
        #### overall ####

        ## baseline hazard from the Cox model
        baseline_Overall <- predictCox(object$models[["OverallSurvival"]], centered = FALSE,
                                       times = eventTimes, newdata = NULL,
                                       type = c("hazard","cumhazard"), 
                                       keep.strata = TRUE, keep.times = TRUE,
                                       se = FALSE, keep.infoVar = TRUE)
        ls.infoVar[["OverallSurvival"]] <- baseline_Cause$infoVar
        
        ## linear predictor for the new observations
        eXb_Overall <- cbind(exp(coxLP(object$models[["OverallSurvival"]], data = newdata, center = FALSE)))
        ## strata for the new observations
        strata_Overall <- coxStrata(object$models[["OverallSurvival"]],
                                    data = newdata,
                                    sterms = ls.infoVar[["OverallSurvival"]]$sterms,
                                    strata.vars = ls.infoVar[["OverallSurvival"]]$strata.vars,
                                    levels = levels(baseline_Overall$strata),
                                    strata.levels = ls.infoVar[["OverallSurvival"]]$strata.levels)
        
        
        #### store
        ls.hazard <- list(matrix(baseline_Cause$hazard, byrow = FALSE, nrow = nEventTimes),# for predictCIF_cpp
                          matrix(baseline_Overall$hazard, byrow = FALSE, nrow = nEventTimes)) # for predictCIF_cpp
        ls.cumhazard <- list(matrix(baseline_Cause$cumhazard, byrow = FALSE, nrow = nEventTimes), # for calcSeCSC
                             matrix(baseline_Overall$cumhazard, byrow = FALSE, nrow = nEventTimes))
        M.eXb <- cbind(eXb_Cause, eXb_Overall)
        M.strata <- cbind(as.numeric(strata_Cause)-1,
                          as.numeric(strata_Overall)-1)
        M.etimes.max <- cbind(baseline_Cause$lastEventTime[M.strata[,1]+1]) # last time by strata
    }

    CIF <- predictCIF_cpp(hazard = ls.hazard, 
                          cumhazard = ls.cumhazard, 
                          eXb = M.eXb, 
                          strata = M.strata,
                          newtimes = sort(times), 
                          etimes = eventTimes, 
                          etimeMax = apply(M.etimes.max,1,min), 
                          t0 = landmark,
                          nEventTimes = nEventTimes,
                          nNewTimes = length(times), 
                          nData = new.n,
                          cause = index.cause - 1, 
                          nCause = nCause,
                          survtype = (surv.type=="survival"),
                          productLimit = product.limit)
    
#### standard error ####
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


        nVar <- unlist(lapply(ls.infoVar,function(m){
            length(m$lpvars)
        }))

        out.seCSC <- calcSeCSC(object,
                               cif = CIF,
                               hazard = ls.hazard,
                               cumhazard = ls.cumhazard,
                               object.time = eventTimes,
                               object.maxtime = apply(M.etimes.max,1,min), 
                               eXb = M.eXb,
                               new.LPdata = new.LPdata,
                               new.strata = M.strata,                               
                               times = sort(times),
                               ls.infoVar = ls.infoVar,
                               new.n = new.n,
                               cause = which(causes == cause),
                               nCause = nCause,
                               nVar = nVar,
                               surv.type = surv.type,
                               export = c("iid"[(iid+band)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE]),
                               store.iid = store.iid)
    }
    
#### export ####
    out <- list(absRisk = CIF[,ootimes,drop=FALSE]) # reorder prediction times

    if(se+band){
        out$absRisk.se <- out.seCSC$se[,ootimes,drop=FALSE]
    }
    if(iid+band){
        out$absRisk.iid <- out.seCSC$iid[,ootimes,,drop=FALSE]
    }
    if(average.iid){
        out$absRisk.average.iid <- out.seCSC$average.iid[,ootimes,drop=FALSE]
    }
    if(keep.times){out$times <- times}

    all.covars <- unique(unlist(lapply(ls.infoVar, function(iI){
        c(iI$lpvars.original, iI$strata.vars.original)
    })))
    if(keep.newdata==TRUE && length(all.covars)>0){
        out$newdata <- newdata[, all.covars, with = FALSE]
    }
    if(keep.strata==TRUE){
        allStrata <- unique(unlist(lapply(ls.infoVar,"[[","strata.vars.original")))
        if (length(allStrata)>0){
            newdata <- copy(newdata[,allStrata, with = FALSE])
            newdata[, (allStrata) := lapply(allStrata, function(col){paste0(col,"=",.SD[[col]])})]
            out$strata <- newdata[, interaction(.SD, sep = " "), .SDcols = allStrata]
        }
    }
    out <- c(out,
             list(se = se,
                  band = band))
    class(out) <- "predictCSC"
    if(confint){
        out <- stats::confint(out)
    }
    if(band && se==FALSE){
        out["absRisk.se"] <- NULL
    }
    if(band && iid==FALSE){
        out["absRisk.iid"] <- NULL
    }
    
    return(out)
}



