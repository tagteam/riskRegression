### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: jun 17 2018 (19:27) 
##           By: Brice Ozenne
##     Update #: 251
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * calcSeCSC (documentation)
#' @title Standard error of the absolute risk predicted from cause-specific Cox models
#' @description  Standard error of the absolute risk predicted from cause-specific Cox models
#' using a first order von Mises expansion of the absolute risk functional.
#' @name calcSeCSC
#' 
#' @param object The fitted cause specific Cox model
#' @param cif the cumulative incidence function at each prediction time for each individual.
#' @param hazard list containing the baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param cumhazard list containing the cumulative baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param object.time a vector containing all the events regardless to the cause.
#' @param object.maxtime a matrix containing the latest event in the strata of the observation for each cause.
#' @param eXb a matrix containing the exponential of the linear predictor evaluated for the new observations (rows) for each cause (columns)
#' @param new.LPdata a list of design matrices for the new observations for each cause.
#' @param new.strata a matrix containing the strata indicator for each observation and each cause.
#' @param times the time points at which to evaluate the predictions.  
#' @param surv.type see the surv.type argument of \code{\link{CSC}}.
#' @param ls.infoVar A list containing the output of \code{coxVariableName} for each Cox model.
#' @param new.n the number of new observations.
#' @param cause the cause of interest.
#' @param nCause the number of causes.
#' @param nVar the number of variables that form the linear predictor in each Cox model
#' @param export can be "iid" to return the value of the influence function for each observation
#'                      "se" to return the standard error for a given timepoint
#' 
#' @param store.iid the method used to compute the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}. See the details section.
#'
#' @details Can also return the empirical influence function of the functionals cumulative hazard or survival
#' or the sum over the observations of the empirical influence function.
#'
#' \code{store.iid="full"} compute the influence function for each observation at each time in the argument \code{times}
#' before computing the standard error / influence functions.
#' \code{store.iid="minimal"} recompute for each subject specific prediction the influence function for the baseline hazard.
#' This avoid to store all the influence functions but may lead to repeated evaluation of the influence function.
#' This solution is therefore efficient more efficient in memory usage but may not be in term of computation time.

## * calcSeCSC (code)
#' @rdname calcSeCSC
calcSeCSC <- function(object, cif, hazard, cumhazard, object.time, object.maxtime,
                      eXb, new.LPdata, new.strata, times, surv.type, ls.infoVar,
                      new.n, cause, nCause, nVar, export, store.iid){

    status <- strata.num <- NULL ## [:CRANcheck:] data.table
                                        # {{{ influence function for each Cox model
    if(is.null(object$iid)){
        object$iid <- list()
        for(iModel in 1:nCause){ # iModel <- 1
            object$iid[[iModel]] <- iidCox(object$models[[iModel]], tau.hazard = object.time, store.iid = store.iid)
        }
    }else{
        store.iid <- object$iid[[1]]$store.iid
        for(iModel in 1:nCause){
            object$iid[[iModel]] <- selectJump(object$iid[[iModel]], times = object.time,
                                               type = c("hazard","cumhazard"))
        }
    }
                                        # }}}
                                        # {{{ prepare arguments
    nEtimes <- length(object.time)
    object.n <- NROW(object$iid[[1]]$IFbeta)
        
    out <- list()
    if("se" %in% export){  
        out$se <- matrix(NA, nrow = new.n, ncol = length(times))
    }    
    if("iid" %in% export){
        out$iid <- array(NA, dim = c(new.n, length(times), object.n))
    }
    if("average.iid" %in% export){
        out$average.iid <- matrix(0, nrow = object.n, ncol = length(times))
    }
                                        # }}}
 
  
    if(store.iid == "minimal"){
        ## {{{ method = "minimal"
        object.modelFrame <- coxModelFrame(object$models[[cause]])
        object.modelFrame[, c("strata.num") := as.numeric(.SD$strata)-1]
        
        ## recover strata in the training dataset
        M.object.strata.num <- do.call(cbind,lapply(1:nCause, function(iCause){
            as.numeric(coxStrata(object$models[[iCause]], data = NULL,
                                 strata.vars = ls.infoVar[[iCause]]$stratavars))-1
        }))
        object.Ustrata <- do.call(paste0, as.data.frame(M.object.strata.num))
        
        ## collapse strata variables into on variable
        level.strata <- unique(new.strata)
        new.Ustrata <- do.call(paste0, as.data.frame(new.strata))
        level.Ustrata <- unique(new.Ustrata)
        nStrata <- length(level.Ustrata)
        
        for(iStrata in 1:nStrata){ # iStrata <- 1

            # {{{ prepare arguments
            nTime <- length(times)
            iStrataTheCause <- new.strata[iStrata,cause]
            indexStrataTheCause <- which(new.Ustrata==level.Ustrata[iStrata])
            
            ls.args <- list()
            ls.args$seqTau <- times
            ls.args$jumpTime <- object.time
            ls.args$jumpTheCause <- object.time %in% object.modelFrame[status==1&strata.num==iStrataTheCause,stop]
            ls.args$IFbeta <- lapply(object$iid, function(x){x$IFbeta})
            ls.args$cif <- cif[indexStrataTheCause,,drop=FALSE]
            ls.args$iS0 <- do.call(cbind,lapply(object$iid, function(x){
                x$calcIFhazard$delta_iS0[[iStrataTheCause+1]]
            }))
            ls.args$newEXb <- eXb[indexStrataTheCause,,drop=FALSE]            
            ls.args$sampleEXb <- exp(do.call(cbind,lapply(object$models, coxLP, data = NULL, center = FALSE)))
            ls.args$sampleTime <- object.modelFrame[["stop"]]

            ls.args$indexJump <- matrix(NA, ncol = nCause, nrow = nEtimes)
            ls.args$indexSample <- matrix(NA, ncol = nCause, nrow = object.n)
            ls.args$Ehazard0 <- vector(mode = "list", length = nCause)
            ls.args$cumEhazard0 <- vector(mode = "list", length = nCause)
            ls.args$cumhazard_iS0 <- vector(mode = "list", length = nCause)
            ls.args$hazard_iS0 <- vector(mode = "list", length = nCause)
            ls.args$X <- vector(mode = "list", length = nCause)
            ls.args$sameStrata <- matrix(NA,nrow = object.n, ncol = nCause)
            ls.args$cumhazard0 <- vector(mode = "list", length = nCause)
            ls.args$hazard0 <- vector(mode = "list", length = nCause)

            for(iterC in 1:nCause){ # iterC <- 1
                iStrataCause <- level.strata[iStrata,iterC]

                ls.args$indexJump[,iterC] <- prodlim::sindex(object$iid[[iterC]]$calcIFhazard$time1[[iStrataCause + 1]],
                                                         eval.times = object.time)
                ls.args$indexSample[,iterC] <- prodlim::sindex(object$iid[[iterC]]$calcIFhazard$time1[[iStrataCause + 1]],
                                                        eval.times = object.modelFrame[["stop"]])
                ls.args$cumEhazard0[[iterC]] <- object$iid[[iterC]]$calcIFhazard$cumElambda0[[iStrataCause + 1]]
                ls.args$Ehazard0[[iterC]] <- object$iid[[iterC]]$calcIFhazard$Elambda0[[iStrataCause + 1]]
                ls.args$cumhazard_iS0[[iterC]] <- c(0,object$iid[[iterC]]$calcIFhazard$cumLambda0_iS0[[iStrataCause + 1]])## add 0 to match prodlim
                ls.args$hazard_iS0[[iterC]] <- c(0,object$iid[[iterC]]$calcIFhazard$lambda0_iS0[[iStrataCause + 1]]) ## add 0 to match prodlim
                ls.args$X[[iterC]] <- new.LPdata[[iterC]][indexStrataTheCause,,drop=FALSE]
                ls.args$sameStrata[,iterC] <- M.object.strata.num[,iterC]==iStrataCause
                ls.args$hazard0[[iterC]] <- hazard[[iterC]][,iStrataCause + 1]
                ls.args$cumhazard0[[iterC]] <- cumhazard[[iterC]][,iStrataCause + 1]
            }
            
            # }}}

            ls.args$theCause <- cause-1
            ls.args$firstJumpTime <- object$iid[[cause]]$etime1.min[iStrataTheCause+1]
            ls.args$lastSampleTime <- object.modelFrame[strata.num==iStrataCause,max(.SD$stop)]
            ls.args$nTau <- nTime
            ls.args$nJump <- nEtimes
            ls.args$nNewObs <- length(indexStrataTheCause)
            ls.args$nSample <- object.n
            ls.args$nCause <- nCause
            ls.args$p <- nVar
            ls.args$survtype <- surv.type=="survival"
            ls.args$exportSE<- ("se" %in% export)
            ls.args$exportIF <- ("iid" %in% export)
            ls.args$exportIFsum <- ("average.iid" %in% export)

            resCpp <- do.call(calcSeCif_cpp, args = ls.args)
            
            if("se" %in% export){
                out$se[indexStrataTheCause,] <- resCpp$se
            }
            if("iid" %in% export){
                out$iid[indexStrataTheCause,,] <- resCpp$iid
            }
            if("average.iid" %in% export){                
                out$average.iid <- out$average.iid + resCpp$iidsum/new.n
            }
            
        }
             
        
        
        # }}}
    }else{
        # {{{ other method
        
        for(iObs in 1:new.n){

            ## compute the individual specific influence function for the cumlative hazard and d(cumulative hazard)
            iStrata <- new.strata[iObs,]        
            iCumHazard <- rep(0, nEtimes)
            iIFhazard1 <- NULL
            iIFcumhazard <- matrix(0, nrow = object.n, ncol = nEtimes)
            for(iCause in 1:nCause){
                X_IFbeta <- tcrossprod(object$iid[[iCause]]$IFbeta,
                                       new.LPdata[[iCause]][iObs,,drop=FALSE])
                if(surv.type == "hazard" || cause != iCause){
                    iCumHazard <- iCumHazard + cumhazard[[iCause]][,iStrata[iCause]+1]*eXb[iObs,iCause]
                                        # Evaluate the influence function for the
                                        # cumulative hazard based on the one of the cumulative baseline hazard
                    if(nVar[iCause] == 0){
                        iIFcumhazard <- iIFcumhazard + object$iid[[iCause]]$IFcumhazard[[iStrata[iCause]+1]]
                    }else{
                        iIFcumhazard <- iIFcumhazard + eXb[iObs,iCause]*(object$iid[[iCause]]$IFcumhazard[[iStrata[iCause]+1]] + crossprod(t(X_IFbeta),cumhazard[[iCause]][,iStrata[iCause]+1]))
                    }
                }
                if(cause == iCause){
                    iHazard1 <- hazard[[cause]][,iStrata[cause]+1]*eXb[iObs,iCause]
                                        # Evaluate the influence function for the
                                        # hazard based on the one of the baseline hazard
                    if(nVar[iCause] == 0){
                        iIFhazard1 <- object$iid[[iCause]]$IFhazard[[iStrata[iCause]+1]]
                    } else{
                        iIFhazard1 <- eXb[iObs,iCause]*(object$iid[[iCause]]$IFhazard[[iStrata[iCause]+1]] + crossprod(t(X_IFbeta),hazard[[iCause]][,iStrata[iCause]+1]))
                    }
                }
            }
        
            
            ## set to s-
            if(nEtimes>1){
                iIFcumhazard <- cbind(0,iIFcumhazard[,1:(nEtimes-1),drop=FALSE])
                iCumHazard <- c(0,iCumHazard[1:(nEtimes-1),drop=FALSE])
            }else{
                iIFcumhazard[] <- 0
                iCumHazard <- 0
            }

            ## influence function for CIF
            IF_tempo <- rowCumSum(rowMultiply_cpp(iIFhazard1 - rowMultiply_cpp(iIFcumhazard, scale = iHazard1),
                                                  scale = exp(-iCumHazard)))    
            IF_tempo <- cbind(0,IF_tempo)[,prodlim::sindex(object.time, eval.times = times)+1,drop=FALSE]
            if(any(times > object.maxtime[iObs])){ # add NA after the last event in the strata
                IF_tempo[,times > object.maxtime[iObs]] <- NA
            
            }

            ## export quantities of interest           
            if("se" %in% export){
                out$se[iObs,] <- sqrt(colSums(IF_tempo^2))
            }
            if("iid" %in% export){
                out$iid[iObs,,] <- t(IF_tempo)
            }
            if("average.iid" %in% export){
                out$average.iid <- out$average.iid + IF_tempo/new.n
            }
            
    }
        # }}}
    }
    
    return(out)
}


#----------------------------------------------------------------------
### calcSeCSC.R ends here
