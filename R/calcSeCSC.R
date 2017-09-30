### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: Sep 30 2017 (18:14) 
##           By: Thomas Alexander Gerds
##     Update #: 167
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
#' @description  Standard error of the absolute risk predicted from cause-specific Cox models
#' using a first order von Mises expansion of the absolute risk functional.
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
#' @param new.n the number of new observations.
#' @param cause the cause of interest.
#' @param nCause the number of causes.
#' @param nVar the number of variables that form the linear predictor in each Cox model
#' @param log.transform Should the variance/influence function be computed on the log(-log) scale
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
#' 
calcSeCSC <- function(object, cif, hazard, cumhazard, object.time, object.maxtime,
                      eXb, new.LPdata, new.strata, times, surv.type,
                      new.n, cause, nCause, nVar, log.transform, export, store.iid){

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

    if(all(c("iid","average.iid") %in% export)){
        stop("Cannot export both average.iid and iid \n")
    }

    ## extract event times
    design <- coxDesign(object$models[[cause]])
    if("strata" %in% names(design) == FALSE){
        design$strata <- 1
    }
    
    out <- list()
    if("se" %in% export){  
        out$se <- matrix(NA, nrow = new.n, ncol = length(times))
    }    
    if("iid" %in% export){
        out$iid <- array(NA, dim = c(new.n, length(times), object.n))
    }
    if("average.iid" %in% export){
        out$iid <- matrix(0, nrow = object.n, ncol = length(times))
    }
    # }}}
 
  
    if(store.iid == "minimal"){
        # {{{ method = "minimal"        
        sample.time <- design[["stop"]]
        ## form all strata
        unique.strata <- unique(new.strata)
        vec.uniqueStrata <- as.numeric(as.factor(apply(unique.strata,1,paste,collapse = ".")))
        nStrata <- length(vec.uniqueStrata)
   
        
        for(iStrata in 1:nStrata){ # iStrata <- 1
            indexStrata <- which(apply(new.strata,1,function(x){
                all(x == unique.strata[iStrata,])
            }))

            # {{{ prepare arguments
            nTime <- length(times)
            iStrataTheCause <- unique.strata[iStrata,1]            
            iStrataTheCause2 <- iStrataTheCause+1
              
            ls.args <- list()
            ls.args$seqTau <- times
            ls.args$jumpTime <- object.time
            ls.args$jumpTheCause <- object.time %in% sample.time[(design[["status"]]==1)*(design[["strata"]]==iStrataTheCause2)==1]
            ls.args$indexJump <- matrix(NA, ncol = nCause, nrow = nEtimes)
            ls.args$indexSample <- matrix(NA, ncol = nCause, nrow = object.n)
            ls.args$IFbeta <- lapply(object$iid, function(x){x$IFbeta})
            ls.args$cif <- cif[indexStrata,,drop=FALSE]
            ls.args$Ehazard0 <- vector(mode = "list", length = nCause)
            ls.args$cumEhazard0 <- vector(mode = "list", length = nCause)
            ls.args$iS0 <- do.call(cbind,lapply(object$iid, function(x){
                x$calcIFhazard$delta_iS0[[iStrataTheCause2]]
            }))
            ls.args$cumhazard_iS0 <- vector(mode = "list", length = nCause)
            ls.args$hazard_iS0 <- vector(mode = "list", length = nCause)
            ls.args$newEXb <- eXb[indexStrata,,drop=FALSE]            
            ls.args$sampleEXb <- exp(do.call(cbind,lapply(object$models, coxLP, data = NULL, center = FALSE)))
            ls.args$X <- vector(mode = "list", length = nCause)
            ls.args$sameStrata <- matrix(NA,nrow = object.n, ncol = nCause)
            ls.args$sampleTime <- sample.time
            ls.args$cumhazard0 <- vector(mode = "list", length = nCause)
            ls.args$hazard0 <- vector(mode = "list", length = nCause)
 
            for(iterC in 1:nCause){ # iterC <- 1
                iStrataCause <- unique.strata[iStrata,iterC]
                iStrataCause2 <- iStrataCause + 1
                
                infoVar <- coxVariableName(object$models[[iterC]])
                object.strata <- as.numeric(coxStrata(object$models[[iterC]], data = NULL,
                                                      strata.vars = infoVar$strata.vars))-1
                ls.args$indexJump[,iterC] <- prodlim::sindex(object$iid[[iterC]]$calcIFhazard$time1[[iStrataCause2]],
                                                         eval.times = object.time)
                ls.args$indexSample[,iterC] <- prodlim::sindex(object$iid[[iterC]]$calcIFhazard$time1[[iStrataCause2]],
                                                        eval.times = sample.time)
                ls.args$cumEhazard0[[iterC]] <- object$iid[[iterC]]$calcIFhazard$cumElambda0[[iStrataCause2]]
                ls.args$Ehazard0[[iterC]] <- object$iid[[iterC]]$calcIFhazard$Elambda0[[iStrataCause2]]
                ls.args$cumhazard_iS0[[iterC]] <- c(0,object$iid[[iterC]]$calcIFhazard$cumLambda0_iS0[[iStrataCause2]])## add 0 to match prodlim
                ls.args$hazard_iS0[[iterC]] <- c(0,object$iid[[iterC]]$calcIFhazard$lambda0_iS0[[iStrataCause2]]) ## add 0 to match prodlim
                ls.args$X[[iterC]] <- new.LPdata[[iterC]][indexStrata,,drop=FALSE]
                ls.args$sameStrata[,iterC] <- object.strata==iStrataCause
                ls.args$hazard0[[iterC]] <- hazard[[iterC]][,iStrataCause2]
                ls.args$cumhazard0[[iterC]] <- cumhazard[[iterC]][,iStrataCause2]
            }
            
            # }}}

            ls.args$theCause <- cause-1
            ls.args$firstJumpTime <- object$iid[[cause]]$etime1.min[iStrataTheCause2]
            ls.args$lastSampleTime <- max(design[design$strata==iStrataCause2,"stop"])
            ls.args$nTau <- nTime
            ls.args$nJump <- nEtimes
            ls.args$nNewObs <- length(indexStrata)
            ls.args$nSample <- object.n
            ls.args$nCause <- nCause
            ls.args$p <- nVar
            ls.args$survtype <- surv.type=="survival"
            ls.args$exportSE<- ("se" %in% export)
            ls.args$exportIF <- ("iid" %in% export)
            ls.args$exportIFsum <- ("average.iid" %in% export)
            ls.args$logTransform <- log.transform

            resCpp <- do.call(calcSeCif_cpp, args = ls.args)
            
            if("se" %in% export){
                out$se[indexStrata,] <- resCpp$se
            }
            if("iid" %in% export){
                out$iid[indexStrata,,] <- resCpp$iid
            }
            if("average.iid" %in% export){
                out$iid <- out$iid + resCpp$iidsum/new.n
            }
            
        }
             
        
        
        # }}}
    }else{
        # {{{ other method
        ## form all strata
        unique.strata <- unique(design$strata)-1
        nStrata <- length(unique.strata)

        ##
        first.event <- sapply(1:nStrata, function(strat){ # strat <- 1
            min(design[(design$status==1)*(design$strata==strat)==1,"stop"])
        })
         
        for(iObs in 1:new.n){
            iStrata <- new.strata[iObs,]        
            iCumHazard <- rep(0, nEtimes)
            iIFhazard1 <- NULL
            iIFcumhazard <- matrix(0, nrow = object.n, ncol = nEtimes)
            
        
        for(iCause in 1:nCause){
            X_IFbeta <- object$iid[[iCause]]$IFbeta %*% t(new.LPdata[[iCause]][iObs,,drop=FALSE])

            if(surv.type == "hazard" || cause != iCause){
                iCumHazard <- iCumHazard + cumhazard[[iCause]][,iStrata[iCause]+1]*eXb[iObs,iCause]
                iIFcumhazard <- iIFcumhazard + IFlambda2hazard(eXb = eXb[iObs,iCause],
                                                               lambda0 = cumhazard[[iCause]][,iStrata[iCause]+1],
                                                               X_IFbeta = X_IFbeta,
                                                               IFlambda0 = object$iid[[iCause]]$IFcumhazard[[iStrata[iCause]+1]],
                                                               nVar = nVar[iCause])
            }
            
            if(cause == iCause){
                iHazard1 <- hazard[[cause]][,iStrata[cause]+1]*eXb[iObs,iCause]

                iIFhazard1 <- IFlambda2hazard(eXb = eXb[iObs,iCause],
                                              lambda0 = hazard[[iCause]][,iStrata[iCause]+1],
                                              X_IFbeta = X_IFbeta,
                                              IFlambda0 = object$iid[[iCause]]$IFhazard[[iStrata[iCause]+1]],
                                              nVar = nVar[iCause])
              
            }
         
        }
        
        # set to s-
        if(nEtimes>1){
            iIFcumhazard <- cbind(0,iIFcumhazard[,1:(nEtimes-1),drop=FALSE])
            iCumHazard <- c(0,iCumHazard[1:(nEtimes-1),drop=FALSE])
        }else{
            iIFcumhazard[] <- 0
            iCumHazard <- 0
        }

            IF_tempo <- rowCumSum(rowMultiply_cpp(iIFhazard1 - rowMultiply_cpp(iIFcumhazard, scale = iHazard1),
                                                  scale = exp(-iCumHazard)))    
            IF_tempo <- cbind(0,IF_tempo)[,prodlim::sindex(object.time, eval.times = times)+1,drop=FALSE]
            if(any(times > object.maxtime[iObs])){ # add NA after the last event in the strata
                IF_tempo[,times > object.maxtime[iObs]] <- NA
            
        }

            if(log.transform){
                IF_tempo <- rowScale_cpp(IF_tempo, scale = cif[iObs,,drop=FALSE]*log(cif[iObs,,drop=FALSE]))
                if(any(times<first.event[new.strata[iObs,cause]+1])){ # any(times<[iObs.strata])
                    IF_tempo[,times<first.event[new.strata[iObs,cause]+1]] <- 0
                }
            }
        
        if("se" %in% export){
            out$se[iObs,] <- sqrt(apply(IF_tempo^2,2,sum))
        }
        if("iid" %in% export){
            out$iid[iObs,,] <- t(IF_tempo)
        }

    }
        # }}}
    }   
    return(out)
}


#----------------------------------------------------------------------
### calcSeCSC.R ends here
