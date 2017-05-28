### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: maj 27 2017 (21:23) 
##           By: Brice Ozenne
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Standard error of the absolute risk predicted from cause-specific Cox models
#' @rdname calcSeCSC
#'
#' @description  Standard error of the absolute risk predicted from cause-specific Cox models.
#'
#' @param cif the cumulative incidence function at each prediction time for each individual.
#' @param hazard list containing the baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param cumhazard list containing the cumulative baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param object.time a vector containing all the events regardless to the cause.
#' @param object.maxtime a matrix containing the latest event in the strata of the observation for each cause.
#' @param iid the value of the influence function for each cause 
#' @param eXb a matrix containing the exponential of the linear predictor evaluated for the new observations (rows) for each cause (columns)
#' @param survtype see the survtype argument of \code{\link{CSC}}.
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
calcSeCSC <- function(cif, hazard, cumhazard, object.time, object.maxtime, iid,
                  eXb, new.LPdata, new.strata, times, survtype,
                  new.n, cause, nCause, nVar, logTransform, export){

    out <- list()
    nEtimes <- length(object.time)
    object.n <- NROW(iid[[1]]$IFbeta)
    if("se" %in% export){  
        out$se <- matrix(NA, nrow = new.n, ncol = length(times))
    }
    if("iid" %in% export){
        out$iid <- array(NA, dim = c(new.n, length(times), object.n))
    }
    
    for(iObs in 1:new.n){
        iStrata <- new.strata[iObs,]        
        iCumHazard <- rep(0, nEtimes)
        iIFhazard1 <- NULL
        iIFcumhazard <- matrix(0, nrow = object.n, ncol = nEtimes)        
        for(iCause in 1:nCause){
            X_IFbeta <- iid[[iCause]]$IFbeta %*% t(new.LPdata[[iCause]][iObs,,drop=FALSE])

            if(survtype == "hazard" || cause != iCause){
                iCumHazard <- iCumHazard + cumhazard[[iCause]][,iStrata[iCause]+1]*eXb[iObs,iCause]
                
                iIFcumhazard <- iIFcumhazard + IFlambda2hazard(eXb = eXb[iObs,iCause],
                                                               lambda0 = cumhazard[[iCause]][,iStrata[iCause]+1],
                                                               X_IFbeta = X_IFbeta,
                                                               IFlambda0 = iid[[iCause]]$IFcumhazard[[iStrata[iCause]+1]],
                                                               nVar = nVar[iCause])
            }
            
            if(cause == iCause){
                iHazard1 <- hazard[[cause]][,iStrata[cause]+1]*eXb[iObs,iCause]
              
                iIFhazard1 <- IFlambda2hazard(eXb = eXb[iObs,iCause],
                                              lambda0 = hazard[[iCause]][,iStrata[iCause]+1],
                                              X_IFbeta = X_IFbeta,
                                              IFlambda0 = iid[[iCause]]$IFhazard[[iStrata[iCause]+1]],
                                              nVar = nVar[iCause])
              
            }
          
        }

        # set to s-
        iIFcumhazard <- cbind(0,iIFcumhazard[,1:(nEtimes-1),drop=FALSE])
        iCumHazard <- c(0,iCumHazard[1:(nEtimes-1),drop=FALSE])

        IF_tempo <- rowCumSum(rowMultiply_cpp(iIFhazard1 - rowMultiply_cpp(iIFcumhazard, scale = iHazard1),
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


#----------------------------------------------------------------------
### calcSeCSC.R ends here
