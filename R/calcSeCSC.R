### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: Jun 25 2020 (14:55) 
##           By: Thomas Alexander Gerds
##     Update #: 956
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
#' @param survival list containing the (all cause) survival in a matrix form. Columns correspond to event times.
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
#' @param diag [logical] when \code{FALSE} the absolute risk/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time. 
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
calcSeCSC <- function(object, cif, hazard, cumhazard, survival, object.time, object.maxtime,
                      eXb, new.LPdata, new.strata, times, surv.type, ls.infoVar,
                      new.n, cause, nCause, nVar, export, store.iid, diag){

    status <- strata.num <- NULL ## [:CRANcheck:] data.table
                                        # {{{ influence function for each Cox model
    if(is.iidCox(object)){
        for(iModel in 1:nCause){ # iModel <- 1
            if(!identical(object$models[[iModel]]$iid$time,object.time)){
                object$models[[iModel]]$iid <- selectJump(object$models[[iModel]]$iid, times = object.time,
                                                          type = c("hazard","cumhazard"))
            }
        }
    }else{
        object <- iidCox(object, tau.hazard = object.time,
                         store.iid = store.iid, return.object = TRUE)
    }
    
    store.iid <- object$models[[1]]$iid$store.iid
                                        # }}}

                                        # {{{ prepare arguments
    nEtimes <- length(object.time)
    object.n <- NROW(object$models[[1]]$iid$IFbeta)
      
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
            ls.args$IFbeta <- lapply(object$models, function(x){x$iid$IFbeta})
            ls.args$iS0 <- do.call(cbind,lapply(object$models, function(x){
                x$iid$calcIFhazard$delta_iS0[[iStrataTheCause+1]]
            }))
            ls.args$newEXb <- eXb[indexStrataTheCause,,drop=FALSE]            
            ls.args$sampleEXb <- exp(do.call(cbind,lapply(object$models, coxLP, data = NULL, center = FALSE)))
            ls.args$sampleTime <- object.modelFrame[["stop"]]
            ls.args$survival <- survival[indexStrataTheCause,,drop=FALSE]
            
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

                ls.args$indexJump[,iterC] <- prodlim::sindex(object$models[[iterC]]$iid$calcIFhazard$time1[[iStrataCause + 1]],
                                                         eval.times = object.time)
                ls.args$indexSample[,iterC] <- prodlim::sindex(object$models[[iterC]]$iid$calcIFhazard$time1[[iStrataCause + 1]],
                                                        eval.times = object.modelFrame[["stop"]])
                ls.args$cumEhazard0[[iterC]] <- object$models[[iterC]]$iid$calcIFhazard$cumElambda0[[iStrataCause + 1]]
                ls.args$Ehazard0[[iterC]] <- object$models[[iterC]]$iid$calcIFhazard$Elambda0[[iStrataCause + 1]]
                ls.args$cumhazard_iS0[[iterC]] <- c(0,object$models[[iterC]]$iid$calcIFhazard$cumLambda0_iS0[[iStrataCause + 1]])## add 0 to match prodlim
                ls.args$hazard_iS0[[iterC]] <- c(0,object$models[[iterC]]$iid$calcIFhazard$lambda0_iS0[[iStrataCause + 1]]) ## add 0 to match prodlim
                ls.args$X[[iterC]] <- new.LPdata[[iterC]][indexStrataTheCause,,drop=FALSE]
                ls.args$sameStrata[,iterC] <- M.object.strata.num[,iterC]==iStrataCause
                ls.args$hazard0[[iterC]] <- hazard[[iterC]][,iStrataCause + 1]
                ls.args$cumhazard0[[iterC]] <- cumhazard[[iterC]][,iStrataCause + 1]
            }
            
            # }}}

            ls.args$theCause <- cause-1
            ls.args$firstJumpTime <- object$models[[cause]]$iid$etime1.min[iStrataTheCause+1]
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

        for(iCause in 1:nCause){ ## remove attributes to have a list of matrices
            attr(new.LPdata[[iCause]],"assign") <- NULL
            attr(new.LPdata[[iCause]],"contrasts") <- NULL
        }

        if("iid" %in% export || "se" %in% export){
            out <- calcSeCif2_cpp(ls_IFbeta = lapply(object$models, function(x){x$iid$IFbeta}),
                                  ls_X = new.LPdata,
                                  ls_cumhazard = cumhazard,
                                  ls_hazard = hazard[[cause]],
                                  survival = survival,
                                  ls_IFcumhazard = lapply(object$models, function(x){x$iid$IFcumhazard}),
                                  ls_IFhazard = object$models[[cause]]$iid$IFhazard,
                                  eXb = eXb,
                                  nJumpTime = nEtimes, JumpMax = object.maxtime,
                                  tau = times, tauIndex = sindex.times, nTau = nTimes,                                  
                                  nObs = object.n,
                                  theCause = (cause-1), nCause = nCause, hazardType = (surv.type=="hazard"), nVar = nVar,
                                  nNewObs = new.n, strata = new.strata,
                                  exportSE = "se" %in% export, exportIF = "iid" %in% export, exportIFsum = "average.iid" %in% export,
                                  diag = diag)
            if("iid" %in% export){
                out$iid <- aperm(out$iid, c(1,3,2))
            }
            
        }else if("average.iid" %in% export){

            ## identify strata index
            new.level.strata <- unique(new.strata)
            new.level.Ustrata <- apply(new.level.strata,1,paste0,collapse="")
            new.n.strata <- length(new.level.Ustrata)

            new.Ustrata <- rowPaste(new.strata)
            new.indexStrata <- lapply(new.level.Ustrata, function(iStrata){
                which(new.Ustrata==iStrata) - 1
            })

            if(is.null(attr(export,"factor"))){
                rm.list <- TRUE
                factor <- vector(mode = "list", length = 1)
            }else{
                rm.list <- FALSE
                factor <- attr(export, "factor")
            }

            all.cause <- switch(surv.type,
                                "hazard" = 1:nCause,
                                "survival" = setdiff(1:nCause,cause))
           
            eMax.obs <- tapply(object.maxtime,new.Ustrata,max)

            test_allCause <- 1:nCause %in% all.cause
            test_theCause <- 1:nCause %in% cause

            n.factor <- length(factor)
            out <- lapply(1:n.factor, function(iF){
                matrix(0, nrow = object.n, ncol = diag + (1-diag)*nTimes)
            })

            ## compute AIF for each strata
            for(iStrata in 1:new.n.strata){ ## iStrata <- 1
                iStrata_causes <- new.level.strata[iStrata,]+1
                iStrata_theCause <- iStrata_causes[cause]
                iIndex_obs <- new.indexStrata[[iStrata]]+1
                iN_obs <- length(iIndex_obs)
                iPrevalence <- iN_obs/new.n

                ## subset times
                if(diag){
                    iTimes <- times[iIndex_obs]
                }else{
                    iTimes <- times
                }
         
                ## subset at jump
                iIndex.jump <- which(hazard[[cause]][,iStrata_theCause]>0)
                iTime.jump <- object.time[iIndex.jump]
                iN.jump <- length(iIndex.jump)

                if(iN.jump==0){next}
                iSindex.times <- prodlim::sindex(jump.times = iTime.jump, eval.times = iTimes)
                iValid.times <- which((iTimes >= min(iTime.jump))*(iTimes <= eMax.obs[iStrata]) == 1)
                iNA.times <- which(iTimes > eMax.obs[iStrata])

                if(diag && length(iNA.times)>0){
                    out <- lapply(out,function(iF){matrix(NA, nrow = nrow(iF), ncol = ncol(iF))})
                    break
                }
                iSindexV.times <- iSindex.times[iValid.times]
                if(diag){ ## only keep individuals and times corresponding to the current strata
                    iIndex_obs <- iIndex_obs[iValid.times]
                    iN_activobs <- length(iIndex_obs)
                }else{
                    iN_activobs <- iN_obs
                }
         
                ## hazard at jump
                ihazard0_cause <- hazard[[cause]][iIndex.jump,iStrata_theCause]
                iIFhazard0_cause <- object$models[[cause]]$iid$IFhazard[[iStrata_theCause]][,iIndex.jump,drop=FALSE]

                ## cumulative hazard just before jump
                iCumhazard0 <- vector(mode = "list", length = nCause)
                iCumhazard0[all.cause] <- lapply(all.cause, function(iC){
                    subsetIndex(cumhazard[[iC]][,iStrata_causes[iC]], index = iIndex.jump-1, default = 0)
                })
                iIFCumhazard0 <- vector(mode = "list", length = nCause)
                iIFCumhazard0[all.cause] <- lapply(all.cause, function(iC){subsetIndex(object$models[[iC]]$iid$IFcumhazard[[iStrata_causes[iC]]],
                                                                                       index = iIndex.jump - 1, default = 0, col = TRUE)})
         
                ## design matrix
                iX <- lapply(1:nCause, function(iC){new.LPdata[[iC]][iIndex_obs,,drop=FALSE]})

                ## linear predictor
                ieXb <- lapply(1:nCause, function(iC){eXb[iIndex_obs, iC]})
                
                ## pre-compute quantities before averaging
                eXb1_S <- colMultiply_cpp(survival[iIndex_obs,iIndex.jump,drop=FALSE], scale = ieXb[[cause]])

                eXb1_S_eXbj <- vector(mode = "list", length = nCause)
                eXb1_S_eXbj[all.cause] <- lapply(all.cause, function(iC){colMultiply_cpp(eXb1_S, ieXb[[iC]])})
                
                if(nVar[cause]>0){
                    eXb1_S_X1 <- lapply(1:nVar[cause], function(iP){ colMultiply_cpp(eXb1_S, scale = iX[[cause]][,iP]) })
                }else{
                    eXb1_S_X1 <- NULL
                }
         
                eXb1_S_Xj_eXbj <- vector(mode = "list", length = nCause)
                for(iCause in all.cause){
                    if(nVar[iCause]>0){
                        eXb1_S_Xj_eXbj[[iCause]] <- lapply(1:nVar[iCause], function(iP){ colMultiply_cpp(eXb1_S_eXbj[[iCause]], scale = iX[[iCause]][,iP]) })
                    }
                }

                for(iFactor in 1:n.factor){
                    ## compute average influence function at each jump time
                    ## <IF>(absRisk) = \int(t) IF_hazard0(cause) E[w * eXb(cause) * S] 
                    ##                       + IF_beta(cause) * hazard0(cause) * E[w * eXb(cause) * X(cause) * S]  
                    ##                       - \sum_j hazard0(cause) * E[w * eXb(cause) * S * eXb(j)] * IF_cumhazard0(j)
                    ##                       - \sum_j hazard0(cause) * E[w * eXb(cause) * X(j) * S * eXb(j)] * cumhazard0(j) * IF_beta(j)

                    if(is.null(factor[[iFactor]])){
                        test.duplicated <- FALSE
                        if(diag){
                            iiN.jump <- iN.jump
                            iiIndex.jump <- 1:iN.jump

                            iMfactor <- do.call(rbind,lapply(iSindexV.times, function(i){c(rep(1,i),rep(0,iN.jump-i))}))
                            iVN_time <- colSums(iMfactor)
                        }else{
                            iiN.jump <- iN.jump
                            iiIndex.jump <- 1:iN.jump

                            iMfactor <- matrix(1, nrow = iN_activobs, ncol = iN.jump)
                            iVN_time <- rep(iN_activobs, iiN.jump)
                        }

                        n.factor2 <- 1
                        
                    }else{
                        if(diag){
                            test.duplicated <- FALSE
                            iiN.jump <- iN.jump
                            iiIndex.jump <- 1:iN.jump

                            iMfactor <- do.call(rbind,lapply(iSindexV.times, function(i){c(rep(1,i),rep(0,iN.jump-i))}))
                            iVN_time <- colSums(iMfactor)
                            iMfactor <- colMultiply_cpp(iMfactor, factor[[iFactor]][iIndex_obs,1])
                            
                            n.factor2 <- 1
                        }else{
                            if(NCOL(factor[[iFactor]]) == 1){ ## same weights at all times
                                iiN.jump <- iN.jump
                                iiIndex.jump <- 1:iN.jump

                                iMfactor <- matrix(factor[[iFactor]][iIndex_obs,1], nrow = iN_activobs, ncol = iN.jump, byrow = FALSE)
                                iVN_time <- rep(iN_activobs, iiN.jump)
                            
                                n.factor2 <- 1
                            }else{ ## weights varying over time
                                n.factor2 <- length(times)
                            }
                        }
                    }


                    for(iFactor2 in 1:n.factor2){ ## iFactor2 <- 1

                        if(n.factor2>1){
                            if((iTimes[iFactor2] < min(iTime.jump)) ||(iFactor2 %in% iNA.times)){ ## 0 or NA, skip (NA are taken care after)
                                next
                            }
                            iiN.jump <- iSindexV.times[iFactor2-sum(iTimes < min(iTime.jump))]
                            iiIndex.jump <- 1:iiN.jump
                            iMfactor <- matrix(factor[[iFactor]][iIndex_obs,iFactor2], nrow = iN_activobs, ncol = iiN.jump, byrow = FALSE)
                            iVN_time <- rep(iN_activobs, iiN.jump)
                        }

                        iAIF <- calcAICcif_R(hazard0_cause = ihazard0_cause,
                                             cumhazard0 = iCumhazard0,
                                             IFhazard0_cause = iIFhazard0_cause,
                                             IFcumhazard0 = iIFCumhazard0,
                                             IFbeta = lapply(object$models,function(iM){iM$iid$IFbeta}),                                         
                                             eXb1_S = eXb1_S,
                                             eXb1_S_eXbj = eXb1_S_eXbj,
                                             eXb1_S_X1 = eXb1_S_X1,
                                             eXb1_S_Xj_eXbj = eXb1_S_Xj_eXbj, 
                                             weight = iVN_time, factor = iMfactor,
                                             nJump = iiN.jump, subsetJump = iiIndex.jump,
                                             nCause = nCause, test_allCause = test_allCause, test_theCause = test_theCause,
                                             nVar = nVar)

             
                        ## export
                        if(diag==TRUE){
                            out[[iFactor]][,1] <- out[[iFactor]][,1] + rowSums(iAIF[,iSindexV.times,drop=FALSE])/iN_obs * iPrevalence
                        }else if(n.factor2>1){
                            out[[iFactor]][,iFactor2] <- out[[iFactor]][,iFactor2] + iAIF[,iiN.jump,drop=FALSE] * iPrevalence
                        }else{
                            out[[iFactor]][,iValid.times] <- out[[iFactor]][,iValid.times] + iAIF[,iSindexV.times] * iPrevalence
                        }
                        
                    }

                    ## NA after last obs
                    if(length(iNA.times)>0){
                        if(diag){
                            out[[iFactor]][iNA.times,] <- NA
                        }else{
                            out[[iFactor]][,iNA.times] <- NA
                        }
                    }

                }
            }

            if(rm.list){
                out <- list(average.iid = out[[1]])
            }else{
                out <- list(average.iid = out)
            }
        
        }
    }

    return(out)
}

## * calcAICcif_R
calcAICcif_R <- function(hazard0_cause,
                         cumhazard0,
                         IFhazard0_cause,
                         IFcumhazard0,
                         IFbeta,                                         
                         eXb1_S,
                         eXb1_S_eXbj,
                         eXb1_S_X1,
                         eXb1_S_Xj_eXbj, 
                         weight, factor,
                         nJump, subsetJump,
                         nCause, test_allCause, test_theCause,
                         nVar){
    
    ## term 1
    iAIF <- rowMultiply_cpp(IFhazard0_cause[,subsetJump,drop=FALSE], scale = colSums(eXb1_S[,subsetJump,drop=FALSE] * factor) / weight)
    for(iCause in 1:nCause){

        ## term 3
        if(test_allCause[iCause]){
            iAIF <- iAIF - rowMultiply_cpp(IFcumhazard0[[iCause]][,subsetJump,drop=FALSE], scale = hazard0_cause[subsetJump] * colSums(eXb1_S_eXbj[[iCause]][,subsetJump,drop=FALSE] * factor) / weight)
        }
             
        if(nVar[iCause]>0){ ## iP <- 1
            E.tempo <- matrix(0, nrow = nVar[iCause], ncol = nJump)

            ## term 2
            if(test_theCause[iCause]){
                ## E.tempo <- E.tempo + do.call(rbind,lapply(eXb1_S_X1, function(iX){colSums(iX[,subsetJump,drop=FALSE] * factor) / weight}))
                for(iP in 1:nVar[iCause]){ ## iP <- 1
                  E.tempo[iP,] <- E.tempo[iP,] + colSums(eXb1_S_X1[[iP]][,subsetJump,drop=FALSE] * factor) / weight
                }
            }
                 
            ## term 4
            if(test_allCause[iCause]){
                ## E.tempo = E.tempo - rowMultiply_cpp(do.call(rbind,lapply(eXb1_S_Xj_eXbj[[iCause]], function(iX){colSums(iX[,subsetJump,drop=FALSE] * factor) / weight})),
                                                    ## cumhazard0[[iCause]][subsetJump])
                for(iP in 1:nVar[iCause]){
                    E.tempo[iP,] <- E.tempo[iP,] - colSums(eXb1_S_Xj_eXbj[[iCause]][[iP]][,subsetJump,drop=FALSE] * factor) * cumhazard0[[iCause]][subsetJump] / weight
                }

            }
            iAIF <- iAIF + IFbeta[[iCause]] %*% rowMultiply_cpp(E.tempo, scale = hazard0_cause[subsetJump])
        }
    }
    
    ## accumulate over time
    return(rowCumSum(iAIF))

}

#----------------------------------------------------------------------
### calcSeCSC.R ends here
