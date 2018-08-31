### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: aug 31 2018 (14:09) 
##           By: Brice Ozenne
##     Update #: 476
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
      
    if(store.iid == "minimal"){
        ## {{{ method = "minimal"
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
        nTimes <- length(times)
        sindex.times <- prodlim::sindex(object.time, eval.times = times)-1 ## i.e. -1 is before the first jump
        if("iid" %in% export || "se" %in% export){            

            for(iCause in 1:nCause){ ## remove attributes to have a list of matrices
                attr(new.LPdata[[iCause]],"assign") <- NULL
                attr(new.LPdata[[iCause]],"contrasts") <- NULL
            }

            out <- calcSeCif2_cpp(ls_IFbeta = lapply(object$iid,"[[","IFbeta"),
                                  ls_X = new.LPdata,
                                  ls_cumhazard = cumhazard,
                                  ls_hazard = hazard[[cause]],
                                  ls_IFcumhazard = lapply(object$iid,"[[","IFcumhazard"),
                                  ls_IFhazard = object$iid[[cause]]$IFhazard,
                                  eXb = eXb,
                                  nJumpTime = nEtimes, JumpMax = object.maxtime,
                                  tau = times, tauIndex = sindex.times, nTau = nTimes,                                  
                                  nObs = object.n,
                                  theCause = (cause-1), nCause = nCause, hazardType = (surv.type=="hazard"), nVar = nVar,
                                  nNewObs = new.n, strata = new.strata,
                                  exportSE = "se" %in% export, exportIF = "iid" %in% export, exportIFsum = "average.iid" %in% export)

            if("iid" %in% export){
                out$iid <- aperm(out$iid, c(1,3,2))
            }
        } else if("average.iid" %in% export){

            new.level.strata <- unique(new.strata)
            new.level.Ustrata <- apply(new.level.strata,1,paste0,collapse="")
            new.n.strata <- length(new.level.Ustrata)
            
            new.Ustrata <- apply(new.strata,1,paste0,collapse="")
            new.indexStrata <- lapply(new.level.Ustrata, function(iStrata){
                which(new.Ustrata==iStrata) - 1
            })

            if(is.null(attr(export,"factor"))){
                rm.list <- TRUE
                factor <- matrix(1, nrow = new.n, ncol = 1)
            }else{
                rm.list <- FALSE
                factor <- attr(export, "factor")
            }
            ## exp(-rowSums(cumhazard[[2]]))
            ## dim(cumhazard[[2]])
            outRcpp <- calcAIFcif_cpp(hazard1 = hazard[[cause]],
                                      ls_cumhazard = cumhazard,
                                      ls_tX = lapply(new.LPdata,t),
                                      eXb = eXb,
                                      ls_IFbeta = lapply(object$iid,"[[","IFbeta"),
                                      ls_IFhazard = object$iid[[cause]]$IFhazard,
                                      ls_IFcumhazard = lapply(object$iid,"[[","IFcumhazard"),
                                      nCause = nCause, theCause = cause-1, hazardType = (surv.type == "hazard"),
                                      nJumpTime = nEtimes, JumpMax = tapply(object.maxtime,new.Ustrata,max),
                                      tau = times, tauIndex = sindex.times, nTau = nTimes,                                  
                                      nObs = object.n, nNewObs = new.n,
                                      levelStrata = new.level.strata, nStrata = new.n.strata, ls_indexStrata = new.indexStrata,
                                      nVar = nVar,
                                      factor = factor)
            
            ## browser()
            if(rm.list){
                out <- list(average.iid = matrix(outRcpp[[1]], nrow = object.n, ncol = nTimes))
            }else{
                out <- list(average.iid <- lapply(outRcpp, function(iMat){matrix(iMat, nrow = object.n, ncol = nTimes)}))
            }

        
        }
                                        # }}}
    }

    return(out)
}

## * old
        ## new.level.strata <- unique(new.strata)
        ## new.level.Ustrata <- apply(new.level.strata,1,paste0,collapse="")
        ## new.n.strata <- length(new.level.Ustrata)
                
        ## new.indexStrata <- lapply(new.level.Ustrata, function(iStrata){
        ## which(new.Ustrata==iStrata) - 1
        ## })
        ## new.prevStrata <- sapply(new.indexStrata, length)/new.n           
 
        ##     if(nCause != 2){
        ##         stop("Argument \'average.iid\' (fast computation) only available when there are exactly two causes.")
        ##     }
        ##     out <- list(average.iid = matrix(NA, nrow = new.n, ncol = nTimes))
        ## if(is.null(attr(export,"factor"))){
        ## rm.list <- TRUE
        ## factor <- list(matrix(1, nrow = new.n, ncol = nTimes))
        ## }else{
        ## rm.list <- FALSE
        ## factor <- attr(export, "factor")
        ## }

        ## browser()
        ## calcAIFcif_cpp(hazard1 = hazard[[cause]],
        ## ls_cumhazard = cumhazard,
        ## ls_X = new.LPdata,
        ## eXb = eXb,
        ## IFhazard1 = object$iid[[cause]]$IFhazard,
        ## ls_IFbeta = lapply(object$iid,"[[","IFbeta"),
        ## ls_IFcumhazard = lapply(object$iid,"[[","IFcumhazard"),
        ## nCause = nCause, theCause = cause - 1, hazardType = (surv.type=="hazard"),
        ## time = times, nTime = nTimes, timeIndex = sindex.times, 
        ## jumpTime = object.time, nJumpTime = nEtimes, maxJumpTime = object.maxtime,
        ## nObs = object.n,
        ## strata = new.strata, nStrata = new.n.strata, prevStrata = new.prevStrata, ls_indexStrata = new.indexStrata,
        ## nVar = nVar,
        ## factor = factor)
        ##     if(nCause != 2){
        ##         stop("Argument \'average.iid\' (fast computation) only available when there are exactly two causes.")
        ##     }
        ##     out <- list(average.iid = matrix(NA, nrow = new.n, ncol = nTimes))
   
        ##     if(nCause != 2){
        ##         stop("Argument \'average.iid\' (fast computation) only available when there are exactly two causes.")
        ##     }
        ##     out <- list(average.iid = matrix(NA, nrow = new.n, ncol = nTimes))
            
        ##     vec.strata <- apply(new.strata,1,paste0,collapse="")
        ##     level.strata <- unique(new.strata)
        ##     level.vec.strata <- apply(level.strata,1,paste0,collapse="")
        ##     n.strata <- length(level.vec.strata)

        ##     n.object.times <- length(object.time)

        ##     for(iStrata in 1:n.strata){ ## iStrata <- 1
        ##         iStrata1 <- level.strata[iStrata,cause] + 1
        ##         iStrata2 <- level.strata[iStrata,-cause] + 1
        ##         iIndex.strata <- which(vec.strata==level.vec.strata[iStrata])
        ##         iN.strata <- length(iIndex.strata)
        ##         iTime.store <- 1
        ##         iPrev.strata <- iN.strata/new.n
                
        ##         int.hazard1 <- 0
        ##         int.b1 <- rep(0, nVar[cause])
        ##         int.b1full <- 0
        ##         int.cumhazard1 <- 0
        ##         int.cumhazard2 <- 0
        ##         int.b2 <- rep(0, nVar[-cause])
        ##         int.b2full <- 0
                
        ##             ## compute integral
        ##             for(iTime in 1:n.object.times){ ## iTime <- 1

        ##                 if(hazard[[cause]][iTime]>0){
        ##                     ## compute expectations at each time
        ##                     E.surv_eXb1 <- 0
        ##                     E.surv_hazard1_X1_1_cumhazard1 <- rep(0, times = nVar[1])
        ##                     E.surv_hazard1_eXb1 <- 0
        ##                     E.surv_hazard1_eXb2 <- 0
        ##                     E.surv_hazard1_X2_cumhazard2 <- rep(0, times = nVar[2])

        ##                     for(iObs in iIndex.strata){ ##  iObs <- 1
        ##                         ## iSurvival <- exp(-sum(sapply(1:nCause, function(iCause){
        ##                             ## cumhazard[[iCause]][iTime]*eXb[iObs,iCause]
        ##                         ## })))
        ##                         iSurvival <- exp(-cumhazard[[1]][iTime]*eXb[iObs,1]-cumhazard[[2]][iTime]*eXb[iObs,2])
        ##                         iHazard1  <- hazard[[cause]][iTime]*eXb[iObs,cause]
        ##                         iCumHazard1  <- cumhazard[[cause]][iTime]*eXb[iObs,cause]
        ##                         iCumHazard2  <- cumhazard[[-cause]][iTime]*eXb[iObs,-cause]
                            
        ##                         E.surv_eXb1 <- E.surv_eXb1 + iSurvival * eXb[iObs,cause]
        ##                         if(nVar[cause]>0){
        ##                             E.surv_hazard1_X1_1_cumhazard1 <- E.surv_hazard1_X1_1_cumhazard1 + iSurvival * iHazard1 * eXb[iObs,cause] * new.LPdata[[cause]][iObs,] * (1 - iCumHazard1)
        ##                         }
        ##                         E.surv_hazard1_eXb1 <- E.surv_hazard1_eXb1 + iSurvival * iHazard1 * eXb[iObs,cause]
        ##                         E.surv_hazard1_eXb2 <- E.surv_hazard1_eXb2 + iSurvival * iHazard1 * eXb[iObs,-cause]
        ##                         if(nVar[-cause]>0){
        ##                             E.surv_hazard1_X2_cumhazard2 <- E.surv_hazard1_X2_cumhazard2 + iSurvival * iHazard1 * new.LPdata[[-cause]][iObs,] * iCumHazard2
        ##                         }
        ##                     }

        ##                     ## update integral
        ##                     int.hazard1 <- int.hazard1 + E.surv_eXb1 * object$iid[[cause]]$IFhazard[[iStrata1]][,iTime]
        ##                     if(nVar[cause]>0){
        ##                         int.b1 <- int.b1 + E.surv_hazard1_X1_1_cumhazard1
        ##                     }
        ##                     int.cumhazard1 <- int.cumhazard1 + E.surv_hazard1_eXb1 * object$iid[[cause]]$IFcumhazard[[iStrata1]][,iTime]
        ##                     int.cumhazard2 <- int.cumhazard2 + E.surv_hazard1_eXb2 * object$iid[[-cause]]$IFcumhazard[[iStrata2]][,iTime]
        ##                     if(nVar[-cause]>0){
        ##                         int.b2 <- int.b2 + E.surv_hazard1_X2_cumhazard2
        ##                     }
        ##                 }
                        
        ##                 ## export
        ##                 if(iTime.store <= nTimes && sindex.times[iTime.store]==iTime){
        ##                     if(nVar[cause]>0){
        ##                         int.b1full <- as.double(object$iid[[cause]]$IFbeta %*% int.b1)
        ##                     }
        ##                     if(nVar[-cause]>0){
        ##                         int.b2full <- as.double(object$iid[[-cause]]$IFbeta %*% int.b2)
        ##                     }
        ##                     ## note: iN.strata corresponds to the 1/n from the expectation of the E. terms
        ##                     int.total <- (int.hazard1 + int.b1full - int.cumhazard1 - int.cumhazard2 - int.b2full) * iPrev.strata / iN.strata
        ##                     out$average.iid[,iTime.store] <- int.total 
        ##                     iTime.store <- iTime.store + 1


        ##                     while(iTime.store <= nTimes && sindex.times[iTime.store]==iTime){
        ##     out$average.iid[,iTime.store] <- int.total
        ##     iTime.store <- iTime.store + 1
        ## }
        ##                 }
                        
        ##             }

                    
        ##     }
                
  
#----------------------------------------------------------------------
### calcSeCSC.R ends here
