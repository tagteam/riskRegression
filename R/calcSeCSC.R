### calcSeCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (21:23) 
## Version: 
## last-updated: sep 20 2019 (18:07) 
##           By: Brice Ozenne
##     Update #: 740
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
calcSeCSC <- function(object, cif, hazard, cumhazard, object.time, object.maxtime,
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
            
        } else if("average.iid" %in% export){

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
                if(diag){
                    warning("Attribute \"factor\" in argument \'average.iid\' is ignored \n")
                }
            }

            all.cause <- switch(surv.type,
                                "hazard" = 1:nCause,
                                "survival" = setdiff(1:nCause,cause))
           
            eMax.obs <- tapply(object.maxtime,new.Ustrata,max)


            out <- warperCalcIID(hazard = hazard, IFhazard = object$models[[cause]]$iid$IFhazard,
                                 cumhazard = cumhazard, IFcumhazard = lapply(object$models,function(iM){iM$iid$IFcumhazard}),
                                 eXb = eXb, X = new.LPdata, IFbeta = lapply(object$models,function(iM){iM$iid$IFbeta}),
                                 cause = cause, all.cause = all.cause, nCause = nCause,
                                 object.n = object.n, new.n = new.n, object.time = object.time, times = times, nTimes = nTimes, 
                                 new.level.strata = new.level.strata, new.n.strata = new.n.strata, new.indexStrata = new.indexStrata,
                                 eMax.obs = eMax.obs, nVar = nVar, factor = factor, diag = diag
                                 )
            out <- list(average.iid = out)

        
        }
                                        # }}}
    }

    return(out)
}

## * warperCalcIID
warperCalcIID <- function(hazard, IFhazard,
                          cumhazard, IFcumhazard,
                          eXb, X, IFbeta,
                          cause, all.cause, nCause,
                          object.n, new.n, object.time, times, nTimes, 
                          new.level.strata, new.n.strata, new.indexStrata,
                          eMax.obs, nVar, diag, factor
                          ){

    test_allCause <- 1:nCause %in% all.cause
    test_theCause <- 1:nCause %in% cause

    n.factor <- length(factor)
    out <- lapply(1:n.factor, function(iF){
        matrix(0, nrow = object.n, ncol = diag + (1-diag)*nTimes)
    })
    
     for(iStrata in 1:new.n.strata){ ## iStrata <- 1
         iStrata_causes <- new.level.strata[iStrata,]+1
         iStrata_theCause <- iStrata_causes[cause]
         iIndex_obs <- new.indexStrata[[iStrata]]+1
         iN_obs <- length(iIndex_obs)
         iPrevalence <- iN_obs/new.n
                
         ## subset at jump
         iIndex.jump <- which(hazard[[cause]][,iStrata_theCause]>0)
         iTime.jump <- object.time[iIndex.jump]
         iN.jump <- length(iIndex.jump)

         if(iN.jump==0){next}
         iSindex.times <- prodlim::sindex(jump.times = iTime.jump, eval.times = times)
         iValid.times <- which((times >= min(iTime.jump))*(times <= eMax.obs[iStrata]) == 1)
         iNA.times <- which(times > eMax.obs[iStrata])

         if(diag && length(iNA.times)>0){
             out <- lapply(out,function(iF){matrix(NA, nrow = nrow(iF), ncol = ncol(iF))})
             break
         }
         iSindexV.times <- iSindex.times[iValid.times]
         if(diag){
             iIndex_obs <- iIndex_obs[iValid.times]
             iN_obs <- length(iIndex_obs)
             iPrevalence <- iN_obs/new.n
         }
         ## hazard at jump
         ihazard0_cause <- hazard[[cause]][iIndex.jump,iStrata_theCause]
         iIFhazard0_cause <- IFhazard[[iStrata_theCause]][,iIndex.jump,drop=FALSE]

         ## cumulative hazard just before jump
         iCumhazard0 <- vector(mode = "list", length = nCause)
         iCumhazard0[all.cause] <- lapply(all.cause, function(iC){
             subsetIndex(cumhazard[[iC]][,iStrata_causes[iC]], index = iIndex.jump-1, default = 0)
         })
         iIFCumhazard0 <- vector(mode = "list", length = nCause)
         iIFCumhazard0[all.cause] <- lapply(all.cause, function(iC){subsetIndex(IFcumhazard[[iC]][[iStrata_causes[iC]]], index = iIndex.jump - 1, default = 0, col = TRUE)})
         
         ## design matrix
         iX <- lapply(1:nCause, function(iC){X[[iC]][iIndex_obs,,drop=FALSE]})

         ## linear predictor
         ieXb <- lapply(1:nCause, function(iC){eXb[iIndex_obs, iC]})

         ## compute survival
         iCumhazard <- lapply(1:nCause, function(iC){
             if(iC %in% all.cause){
                 return(tcrossprod(ieXb[[iC]],iCumhazard0[[iC]]))
             }else{
                 return(matrix(0,iN_obs,iN.jump))
             }
         })
         iSurvival <- exp(-do.call("+",iCumhazard))

         ## pre-compute quantities before averaging
         eXb1_S <- colMultiply_cpp(iSurvival, scale = ieXb[[cause]])

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

             ## ## R version
             if(is.null(factor[[iFactor]])){                 
                 if(diag){
                     iMfactor <- do.call(rbind,lapply(iSindexV.times, function(i){c(rep(1,i),rep(0,iN.jump-i))}))
                     ## iTable <- table(iSindexV.times)
                     iVN <- colSums(iMfactor)
                 }else{
                     iMfactor <- matrix(1, nrow = iN_obs, ncol = iN.jump)
                     iVN <- rep(1,iN_obs)
                 }                 
             }else{
                 iMfactor <- matrix(factor[[iFactor]][,1], nrow = iN_obs, ncol = iN.jump, byrow = FALSE)
                 iVN <- rep(1,iN_obs)
             }

             ## term 1
             iAIF <- rowMultiply_cpp(iIFhazard0_cause, scale = colSums(eXb1_S * iMfactor)/iVN)
             browser()
             for(iCause in 1:nCause){

                 ## term 3
                 if(test_allCause[iCause]){
                     iAIF <- iAIF - rowMultiply_cpp(iIFCumhazard0[[iCause]], scale = ihazard0_cause * colSums(eXb1_S_eXbj[[iCause]] * iMfactor)/iVN)
                 }
             
                 if(nVar[iCause]>0){ ## iP <- 1
                     E.tempo <- rep(0, length = iN.jump)

                     if(test_theCause[iCause]){
                         ## term 2
                         E.tempo = E.tempo + do.call(rbind,lapply(eXb1_S_X1, function(iX){colSums(iX * iMfactor)/iVN}))
                     }
                 
                     ## term 4
                     if(test_allCause[iCause]){
                         E.tempo = E.tempo - rowMultiply_cpp(do.call(rbind,lapply(eXb1_S_Xj_eXbj[[iCause]], function(iX){colSums(iX * iMfactor)/iVN})),
                                                             iCumhazard0[[iCause]])
                     }
                 
                     iAIF <- iAIF + IFbeta[[iCause]] %*% (E.tempo * ihazard0_cause)
                 }
             }
             ## accumulate over time
             ## head(iAIF)
             iAIF <- rowCumSum(iAIF)

             ## export
             if(diag){
                 out[[iFactor]][,1] = out[[iFactor]][,1] + rowMeans(iAIF) * iPrevalence
             }else{
                 out[[iFactor]][,iValid.times] = out[[iFactor]][,iValid.times] + iAIF[,iSindexV.times] * iPrevalence
             }

             ## NA after last obs
             if(length(iNA.times)){
                 if(diag){
                     out[[iFactor]][iNA.times,] <- NA
                 }else{
                     out[[iFactor]][,iNA.times] <- NA
                 }
             }

         }
         ## ## C++ version
         ## xx <- calcAIFcif_cpp(hazard1 = ihazard0_cause,
         ##                      IFhazard1 = iIFhazard0_cause,
         ##                      cumhazard = iCumhazard0,
         ##                      IFcumhazard = iIFCumhazard0,
         ##                      IFbeta = IFbeta,
         ##                      eXb1_S = ieXb,
         ##                      eXb1_S_eXbj = eXb1_S_eXbj,
         ##                      eXb1_S_X1 = eXb1_S_X1,
         ##                      eXb1_S_Xj_eXbj = eXb1_S_Xj_eXbj,
         ##                      nObs = object.n, nJump = length(object.time), nTimes = nTimes,
         ##                      validTimes = iValid.times, NATimes = iNA.times, sindexTimes = iSindex.times, 
         ##                      nCause = nCause, test_allCause = test_allCause, test_theCause = test_theCause,
         ##                      nVar = nVar,
         ##                      factor = list(matrix(1, nrow = new.n, ncol = 1)),
         ##                      diag = FALSE, is_factor2 = FALSE)

     }

     return(out)
}

#----------------------------------------------------------------------
### calcSeCSC.R ends here
