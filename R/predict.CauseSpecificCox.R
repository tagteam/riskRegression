#' @title Predicting absolute risk from cause-specific Cox models
#' @rdname predict.CauseSpecificCox
#' @aliases predict.CauseSpecificCox
#' @aliases predictBig.CauseSpecificCox
#'
#' @description  Apply formula to combine two or more Cox models into absolute risk (cumulative incidence function)
#' 
#' @param object The fitted cause specific Cox model
#' @param newdata A data frame containing the values of the variables
#'     in the right hand side of 'coxph' for each subject.
#' @param times Vector of times at which to return the estimated
#'     hazards/survival
#' @param cause Identifies the cause of interest among the competing
#'     events.
#' @param landmark the starting time for the computation of the cumulative risk
#' @param keep.times Logical. If \code{TRUE} add the evaluation times
#'     to the output.
#' @param keep.newdata Logical. If \code{TRUE} add the value of the covariates used to make the prediction in the output list. 
#' @param keep.strata Logical. If \code{TRUE} add the value of the strata used to make the prediction in the output list. 
#' @param se Logical. If \code{TRUE} add the standard errors to the output.
#' @param band Logical. If \code{TRUE} add the confidence band to the output.
#' @param iid Logical. If \code{TRUE} add the influence function to the output.
#' @param nSim.band the number of simulations used to compute the quantiles
#' for the confidence bands.
#' @param conf.level Level of confidence.
#' @param ci.logTransform Should the confidence intervals/bands be computed on the
#' log(-log) scale and be backtransformed.
#' Otherwise they are computed on the original scale and truncated (if necessary).
#' @param ... not used
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds
#'     tag@@biostat.ku.dk
#' @details Note: for Cox regression models with time varying
#'     covariates it does not make sense to use this function, because
#'     the predicted risk has to be a measurable function of the data
#'     available at the time origin.
#' 
#' @return A list containing:
#' \itemize{
#' \item{absRisk}: (data table) the predictions for each subject (in rows) and each time (in columns).
#' \item{absRisk.se}: (data table) the standard errors of the predictions.
#' \item(absRisk.iid): (array) the value of the influence of each subject used to fit the object (dim 3)
#' for each subject in newdata (dim 1) and each time (dim 2).
#' \item{times}: (vector) the evaluation times.
#' }
#'  
#' @examples 
#' set.seed(5)
#' d <- sampleData(80,outcome="comp")
#' nd <- sampleData(4,outcome="comp")
#' d$time <- round(d$time,1)
#' ttt <- sort(sample(x = unique(d$time), size = 10))
#'
#' # coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X3+X8,data=d, method = "breslow")
#' x= predict(CSC.fit,newdata=nd,times=1:10,cause=1,se=1L)
#' px=print(x)
#' px
#' x2 = predict(CSC.fit,newdata=nd,times=1:10,cause=1,se=1L,
#'            ci.logTransform = TRUE)
#' 
#' predCSC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
#' predCSC.se <- predict(CSC.fit, newdata = d[1:5,], cause = 2, times = ttt,
#'                       se = TRUE,keep.newdata=TRUE)
#' predCSC.iid <- predict(CSC.fit, newdata = d[1:5,],
#'                        cause = 2, times = ttt, iid = TRUE)
#'
#' # predCSC.se$absRisk.se
#' # sqrt(apply(predCSC.iid$absRisk.iid[,1,]^2,1,function(x){sum(x)}))
#' ## strata
#' CSC.fit.s <- CSC(list(Hist(time,event)~ strata(X1)+X2+X9,
#'  Hist(time,event)~ X2+strata(X4)+X8+X7),data=d, method = "breslow")
#' predict(CSC.fit.s,cause=1,times=ttt,se=1L)
#' # cph function
#' CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", fitter = "cph")
#' 
#' predict(CSC.cph, newdata = d, cause = 2, times = ttt)
#' 
#' # landmark analysis
#' T0 <- 1
#' predCSC_afterT0 <- predict(CSC.fit, newdata = d, cause = 2, times = ttt[ttt>T0], landmark = T0)
#' predCSC_afterT0
#'
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
                                     se  = FALSE,
                                     band = FALSE,
                                     iid = FALSE,
                                     nSim.band = 1e4,
                                     ci.logTransform = FALSE,
                                     conf.level=0.95,
                                     ...){
    
    if(object$fitter=="phreg"){newdata$entry <- 0} 
    if(missing(newdata)){newdata <- eval(object$call$data)}
    data.table::setDT(newdata)
    
    survtype <- object$survtype
    if (length(cause) > 1){
        stop(paste0("Can only predict one cause. Provided are: ", 
                    paste(cause, collapse = ", "), sep = ""))
    }
    if (missing(cause)) {
        cause <- object$theCause
    }
    causes <- object$causes
    # cannot use only eventtimes of cause 1 otherwise wrong interpolation in the C++ function
    eTimes <- object$eventTimes
    
    if (any(match(as.character(cause), causes, nomatch = 0)==0L))
        stop(paste0("Requested cause ",as.character(cause)," does not match fitted causes which are:\n ",paste0("- ",causes,collapse="\n")))
    ## stopifnot(match(as.character(cause), causes, nomatch = 0) != 
    ## 0)
    if (survtype == "survival") {
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

    ## Confidence bands
    if(band>0){
        iid.save <- iid
        se.save <- se
        iid <- TRUE
        se <- TRUE
    }
    
    # relevant event times to use  
    eventTimes <- eTimes[which(eTimes <= max(times))] 
    if(length(eventTimes) == 0){eventTimes <- min(times)} # at least the first event

    # order prediction times
    ootimes <- order(order(times))

    # predict cumulative cause specific hazards
    new.n <- NROW(newdata)
    nEventTimes <- length(eventTimes)
    
    if (survtype == "hazard") {
        nCause <- length(causes)
        
        ls.hazard <- vector(mode = "list", length = nCause)
        ls.cumhazard <- vector(mode = "list", length = nCause)
        M.eXb_h <- matrix(NA, nrow = new.n, ncol = nCause)
        M.strata <- matrix(NA, nrow = new.n, ncol = nCause)
        M.etimes.max <- matrix(NA, nrow = new.n, ncol = nCause)
        
        for(iterC in 1:nCause){
            infoVar <- CoxVariableName(object$models[[iterC]])
            
            ## baseline hazard from the Cox model
            causeBaseline <- predictCox(object$models[[iterC]], centered = FALSE,
                                        times = eventTimes, newdata = NULL,
                                        type = c("hazard","cumhazard"), 
                                        keep.strata = TRUE, keep.times = TRUE,
                                        se = FALSE)
            
            ls.hazard[[iterC]] <- matrix(causeBaseline$hazard, byrow = FALSE, nrow = nEventTimes)
            ls.cumhazard[[iterC]] <- matrix(causeBaseline$cumhazard, byrow = FALSE, nrow = nEventTimes) 
            
            ## linear predictor for the new observations
            newdata.eXb <- exp(CoxLP(object$models[[iterC]], data = newdata, center = FALSE))
            
            ## strata for the new observations
            newdata.strata <- CoxStrata(object$models[[iterC]], data = newdata, 
                                        sterms = infoVar$sterms, 
                                        stratavars = infoVar$stratavars, 
                                        levels = levels(causeBaseline$strata), 
                                        stratalevels = infoVar$stratalevels)
            
            M.strata[,iterC] <- as.numeric(newdata.strata)-1
            M.eXb_h[,iterC] <- newdata.eXb
            
            ## last time by strata
            M.etimes.max[,iterC] <- causeBaseline$lastEventTime[M.strata[,iterC]+1]
        }
        
        M.eXb_cumH <- M.eXb_h
        
    }else{
        nCause <- 1
        tdiff <- 0*min(diff(eTimes))/2 # TO MATCH test-CauseSpecificCoxRegresion.R but will not match pec
        #### cause ####
        infoVar_Cause <- CoxVariableName(object$models[[paste("Cause",cause)]])
        
        ## baseline hazard from the Cox model
        causeBaseline <- predictCox(object$models[[paste("Cause",cause)]],
                                    centered = FALSE,
                                    times = eventTimes,
                                    newdata = NULL,
                                    type = c("hazard","cumhazard"),
                                    keep.strata = TRUE,
                                    keep.times = TRUE,
                                    se = FALSE)

        ls.hazard <- list(matrix(causeBaseline$hazard, byrow = FALSE, nrow = nEventTimes))
        
        ## linear predictor for the new observations
        newdata.eXb_cause <- exp(CoxLP(object$models[[paste("Cause",cause)]], data = newdata, center = FALSE))
        M.eXb_h <- cbind(newdata.eXb_cause)
        
        ## strata for the new observations
        newdata.strata <- CoxStrata(object$models[[paste("Cause",cause)]],
                                    data = newdata,
                                    sterms = infoVar_Cause$sterms,
                                    stratavars = infoVar_Cause$stratavars,
                                    levels = levels(causeBaseline$strata),
                                    stratalevels = infoVar_Cause$stratalevels)
        
        M.strata <- cbind(as.numeric(newdata.strata)-1)
        M.etimes.max <- cbind(causeBaseline$lastEventTime[M.strata+1]) # last time by strata
        
        #### overall ####

        ## baseline
        overallBaseline <- predictCox(object$models[["OverallSurvival"]], centered = FALSE,
                                      times = eventTimes-tdiff, newdata = NULL,
                                      type = "cumhazard", 
                                      keep.strata = TRUE, keep.times = TRUE,
                                      se = FALSE)
        ls.cumhazard <- list(matrix(overallBaseline$cumhazard, byrow = FALSE, nrow = nEventTimes))
        
        ## linear predictor for the new observations
        newdata.eXb_All <- exp(CoxLP(object$models[["OverallSurvival"]], data = newdata, center = FALSE))
        M.eXb_cumH <- cbind(newdata.eXb_All)
        
    }

    CIF <- predictCIF_cpp(hazard = ls.hazard, 
                          cumhazard = ls.cumhazard, 
                          eXb_h = M.eXb_h, 
                          eXb_cumH = M.eXb_cumH, 
                          strata = M.strata,
                          newtimes = sort(times), 
                          etimes = eventTimes, 
                          etimeMax = apply(M.etimes.max,1,min), 
                          t0 = landmark,
                          nEventTimes = nEventTimes,
                          nNewTimes = length(times), 
                          nData = new.n,
                          cause = which(causes == cause) - 1, 
                          nCause = nCause)
    
    #### standard error ####
    if(se || iid){
        if(!is.na(landmark)){
            stop("standard error for the conditional survival not implemented \n")
        }

        ## design matrix
        new.LPdata <- list()
        for(iCause in 1:nCause){
            
            infoVar <- CoxVariableName(object$models[[iterC]])
            
            if(length(infoVar$lpvars) > 0){
                new.LPdata[[iCause]] <- model.matrix(object$models[[iCause]], newdata)
            }else{
                new.LPdata[[iCause]] <- matrix(0, ncol = 1, nrow = new.n)
            }  
        }
        
        ## influence function 
        
        if(is.null(object$iid)){
            object$iid <- list()
            for(iModel in 1:nCause){
                object$iid[[iModel]] <- iidCox(object$models[[iModel]], tauHazard = eventTimes)
            }
        }else{
            for(iModel in 1:nCause){
                object$iid[[iModel]] <- selectJump(object$iid[[iModel]], times = eventTimes,
                                                   type = c("hazard","cumhazard"))
            }
        }

        nVar <- unlist(lapply(object$models,function(m){
            length(CoxVariableName(m)$lpvars)
        }))
        
        out.seCSC <- seCSC(cif = CIF,
                           hazard = ls.hazard,
                           cumhazard = ls.cumhazard,
                           object.time = eventTimes,
                           object.maxtime = apply(M.etimes.max,1,min), 
                           iid =  object$iid,
                           eXb_h = M.eXb_h,
                           eXb_cumH = M.eXb_cumH,
                           new.LPdata = new.LPdata,
                           new.strata = M.strata,
                           times = sort(times),
                           new.n = new.n,
                           cause = which(causes == cause),
                           nCause = nCause,
                           nVar = nVar,
                           logTransform = ci.logTransform,
                           export = c("iid"[iid==TRUE],"se"[se==TRUE]))
    }
    
    #### export ####
    out <- list(absRisk = CIF[,ootimes,drop=FALSE]) # reorder prediction times

    if(se){
        out$absRisk.se <- out.seCSC$se[,ootimes,drop=FALSE]
        zval <- qnorm(1-(1-conf.level)/2, 0,1)

        if(ci.logTransform){
            out$absRisk.lower <- exp(-exp(log(-log(out$absRisk)) + zval*out$absRisk.se))
            out$absRisk.upper <- exp(-exp(log(-log(out$absRisk)) - zval*out$absRisk.se))
        }else{            
            # to keep matrix format even when out$absRisk contains only one line
            out$absRisk.lower <- out$absRisk.upper <- matrix(NA, nrow = NROW(out$absRisk.se), ncol = NCOL(out$absRisk.se))
            out$absRisk.lower[] <- apply(out$absRisk - zval*out$absRisk.se,2,pmax,0)
            out$absRisk.upper[] <- apply(out$absRisk + zval*out$absRisk.se,2,pmin,1)
        }
    }
    if(iid){
        out$absRisk.iid <- out.seCSC$iid[,ootimes,,drop=FALSE]
    }
    if(band>0){
        
        out$quantile.band <- confBandCox(iid = out$absRisk.iid,
                                         se = out$absRisk.se,
                                         times = times,                                         
                                         n.sim = nSim.band,
                                         conf.level = conf.level)
            
        if(iid.save==FALSE){
            out$absRisk.iid <- NULL
        }

        quantile95 <- colMultiply_cpp(out$absRisk.se,out$quantile.band)
                
        if(ci.logTransform){
            out$absRisk.lowerBand <- exp(-exp(log(-log(out$absRsik)) + quantile95))
            out$absRisk.upperBand <- exp(-exp(log(-log(out$absRsik)) - quantile95))
        }else{            
            out$absRisk.lowerBand <- matrix(NA, nrow = NROW(out$absRisk.se), ncol = NCOL(out$absRisk.se))
            out$absRisk.lowerBand[] <- apply(out$absRisk - quantile95,2,pmax,0)
            out$absRisk.upperBand <- out$absRisk + quantile95
        }
        
        if(se.save==FALSE){
            out$absRisk.se <- NULL
            out$absRisk.lower <- NULL
            out$absRisk.upper <- NULL
        }
    }    
    if(keep.times){out$times <- times}
    if(keep.newdata==TRUE){
        allVars <- unique(unlist(lapply(object$models, function(m){CoxCovars(m)})))
        if(length(allVars)>0){
            out$newdata <- newdata[, allVars, with = FALSE]
        }
    }
    if(keep.strata==TRUE){
        allStrata <- unique(unlist(lapply(object$models, function(m){CoxVariableName(m)$stratavars.original})))
        if (length(allStrata)>0){
            newdata <- copy(newdata[,allStrata, with = FALSE])
            newdata[, (allStrata) := lapply(allStrata, function(col){paste0(col,"=",.SD[[col]])})]
            out$strata <- newdata[, interaction(.SD, sep = " "), .SDcols = allStrata]
        }
    }
    out$conf.level <- conf.level
    out <- c(out,list(se=se))
    class(out) <- "predictCSC"
    return(out)
}



#' @title Standard error of the absolute risk predicted from cause-specific Cox models
#' @rdname seCSC
#'
#' @description  Standard error of the absolute risk predicted from cause-specific Cox models.
#'
#' @param cif the cumulative incidence function at each prediction time for each individual.
#' @param hazard list containing the baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param cumhazard list containing the cumulative baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param object.time a vector containing all the events regardless to the cause.
#' @param object.maxtime a matrix containing the latest event in the strata of the observation for each cause.
#' @param iid the value of the influence function for each cause 
#' @param eXb_h a matrix containing the exponential of the linear predictor evaluated for the new observations (rows) for each cause (columns)
#' @param eXb_cumH same as before except when considering \code{survtype == "survival"}
#' @param new.LPdata a list of design matrices for the new observations for each cause.
#' @param new.strata a matrix containing the strata indicator for each observation and each cause.
#' @param times the time points at which to evaluate the predictions.  
#' @param new.n the number of new observations.
#' @param cause the cause of interest.
#' @param nCause the number of causes.
#' @param nVar the number of variables that form the linear predictor in each Cox model
#' @param logTransform Should the variance/influence function be computed on the log(-log) scale
#' @param export can be "iid" to return the value of the influence function for each observation
#'                      "se" to return the standard error for a given timepoint
#' 
#' @examples 
#' 
#' set.seed(10)
#' d <- SimCompRisk(2e1)
#' d$time <- round(d$time,1)
#' ttt <- unique(sort(d$time))#sort(sample(x = unique(sort(d$time)), size = 10))
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow")
#' 
#' predCSC <- predict(CSC.fit, newdata = d[1,,drop=FALSE], cause = 2, times = ttt, se = TRUE)
#'
#' 
seCSC <- function(cif, hazard, cumhazard, object.time, object.maxtime, iid,
                  eXb_h, eXb_cumH, new.LPdata, new.strata, times,
                  new.n, cause, nCause, nVar, logTransform, export){

    out <- list()
    nEtimes <- length(object.time)
    object.n <- NROW(iid[[1]]$ICbeta)
    if("se" %in% export){  
        out$se <- matrix(NA, nrow = new.n, ncol = length(times))
    }
    if("iid" %in% export){
        out$iid <- array(NA, dim = c(new.n, length(times), object.n))
    }
    
    for(iObs in 1:new.n){
        iStrata <- new.strata[iObs,]        
        iCumHazard <- rep(0, nEtimes)
        iIChazard1 <- NULL
        iICcumhazard <- matrix(0, nrow = object.n, ncol = nEtimes)    
        for(iCause in 1:nCause){
            #
            iCumHazard <- iCumHazard + cumhazard[[iCause]][,iStrata[iCause]+1]*eXb_cumH[iObs,iCause]
            #
            X_ICbeta <- iid[[iCause]]$ICbeta %*% t(new.LPdata[[iCause]][iObs,,drop=FALSE])

            iICcumhazard <- iICcumhazard + IClambda2hazard(eXb = eXb_cumH[iObs,iCause],
                                                           lambda0 = cumhazard[[iCause]][,iStrata[iCause]+1],
                                                           X_ICbeta = X_ICbeta,
                                                           IClambda0 = iid[[iCause]]$ICcumhazard[[iStrata[iCause]+1]],
                                                           nVar = nVar[iCause])
          
            if(cause == iCause){
                iHazard1 <- hazard[[cause]][,iStrata[cause]+1]*eXb_h[iObs,iCause]
              
                iIChazard1 <- IClambda2hazard(eXb = eXb_h[iObs,iCause],
                                              lambda0 = hazard[[iCause]][,iStrata[iCause]+1],
                                              X_ICbeta = X_ICbeta,
                                              IClambda0 = iid[[iCause]]$IChazard[[iStrata[iCause]+1]],
                                              nVar = nVar[iCause])
              
          }
          
      }

        IF_tempo <- rowCumSum(rowMultiply_cpp(iIChazard1 - rowMultiply_cpp(iICcumhazard, scale = iHazard1),
                                              scale = exp(-iCumHazard)))        
        IF_tempo <- cbind(0,IF_tempo)[,prodlim::sindex(object.time, eval.times = times)+1,drop=FALSE]
        if(any(times > object.maxtime[iObs])){ # add NA after the last event in the strata
             IF_tempo[,times > object.maxtime[iObs]] <- NA
        }

        if(logTransform){
            IF_tempo <- rowScale_cpp(IF_tempo, scale = cif[iObs,,drop=FALSE]*log(cif[iObs,,drop=FALSE])) 
        }
        
        if("se" %in% export){
            out$se[iObs,] <- sqrt(apply(IF_tempo^2,2,sum))
        }
        if("iid" %in% export){
            out$iid[iObs,,] <- t(IF_tempo)
        }

    }
    return(out)
}
