### calcSeCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (11:46) 
## Version: 
## last-updated: May 30 2018 (14:03) 
##           By: Thomas Alexander Gerds
##     Update #: 260
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ calcSeCox

#' Computation of standard errors for predictions
#'
#' Compute the standard error associated to the predictions from Cox regression model
#' using a first order von Mises expansion of the functional (cumulative hazard or survival).
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param times Vector of times at which to return the estimated
#'      hazard/survival.
#' @param nTimes the length of the argument \code{times}. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param Lambda0 the baseline hazard estimate returned by \code{BaseHazStrata_cpp}.
#' @param object.n the number of observations in the dataset used to estimate the object. 
#' @param object.time the time to event of the observations used to estimate the object.
#' @param object.eXb the exponential of the linear predictor relative to the observations used to estimate the object. 
#' @param object.strata the strata index of the observations used to estimate the object.
#' @param nStrata the number of strata.
#' @param new.eXb the linear predictor evaluated for the new observations
#' @param new.LPdata the variables involved in the linear predictor for the new observations
#' @param new.strata the strata indicator for the new observations
#' @param new.survival the survival evaluated for the new observations
#' @param new.cumhazard the cumulative hazard evaluated for the new observations
#' @param nVar the number of variables that form the linear predictor
#' @param log.transform Should the variance/influence function be computed on the log or log(-log) scale
#' @param export can be "iid" to return the value of the influence function for each observation
#'                      "se" to return the standard error for a given timepoint
#'                      
#' @param store.iid Implementation used to estimate the influence function and the standard error.
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
#' @inheritParams  predict.CauseSpecificCox
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
calcSeCox <- function(object, times, nTimes, type, 
                      Lambda0, object.n, object.time, object.eXb, object.strata, nStrata,
                      new.eXb, new.LPdata, new.strata, new.survival, new.cumhazard,
                      nVar, log.transform, export, store.iid){

    
    if("average.iid" %in% export && log.transform){
            stop("Cannot compute average.iid when log.transform equals to TRUE \n")
    }
    
    # {{{ computation of the influence function
    iid.object <- object$iid
    if(is.null(iid.object)){
        iid.object <- iidCox(object, tau.hazard = times, store.iid = store.iid)
    }else{
        store.iid <- iid.object$store.iid
        iid.object <- selectJump(iid.object, times = times, type = type)        
    }
    # }}}

    # {{{ prepare arguments
    n.new <- length(new.eXb)
    new.strata <- as.numeric(new.strata)
  
    if(length(Lambda0$strata)==0){
        Lambda0$strata <- rep(1, length(Lambda0$time))
    }else{
        Lambda0$strata <- as.numeric(Lambda0$strata)    
    }
  
    if("hazard" %in% type){Lambda0$hazard <- lapply(1:nStrata,function(s){Lambda0$hazard[Lambda0$strata==s]})}
    if("cumhazard" %in% type || "survival" %in% type){Lambda0$cumhazard <- lapply(1:nStrata,function(s){Lambda0$cumhazard[Lambda0$strata==s]})}
    
    out <- list()
    if("se" %in% export){
        if("cumhazard" %in% type){out$cumhazard.se <- matrix(NA, nrow = n.new, ncol = nTimes)}
        if("survival" %in% type){out$survival.se <- matrix(NA, nrow = n.new, ncol = nTimes)}
    }
    if("iid" %in% export){
        if("hazard" %in% type){out$hazard.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
        if("cumhazard" %in% type){out$cumhazard.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
        if("survival" %in% type){out$survival.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
    }
    if("average.iid" %in% export){
        if("cumhazard" %in% type){out$cumhazard.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)}
        if("survival" %in% type){out$survival.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)}
    }
    # }}}

    if(store.iid == "minimal"){
        # {{{ method = minimal
        object.strata <- as.numeric(object.strata)
        if("hazard" %in% type){
            stop("store.iid=\"minimal\" cannot be used to extract the influence function of the hazard \n")
        }
        
        for(iStrata in 1:nStrata){ # iStrata <- 1
            indexStrata <- which(new.strata==iStrata)
            if(length(indexStrata)==0){next}

            resCpp <- calcSeHazard_cpp(seqTau = times,
                                       indexTau = prodlim::sindex(iid.object$calcIFhazard$time1[[iStrata]], eval.times = times),
                                       indexJump = prodlim::sindex(iid.object$calcIFhazard$time1[[iStrata]], eval.times = object.time),
                                       IFbeta = iid.object$IFbeta,
                                       cumEhazard0 = iid.object$calcIFhazard$cumElambda0[[iStrata]],
                                       iS0 = iid.object$calcIFhazard$delta_iS0[[iStrata]],
                                       cumhazard_iS0 = c(0,iid.object$calcIFhazard$cumLambda0_iS0[[iStrata]]),
                                       newEXb = new.eXb[indexStrata],
                                       sampleEXb = object.eXb,
                                       X = new.LPdata[indexStrata,,drop=FALSE],
                                       sameStrata = (object.strata==iStrata),
                                       sampleTime = object.time,
                                       cumhazard0 = Lambda0$cumhazard[[iStrata]],
                                       newHazard = new.cumhazard[indexStrata,,drop=FALSE],
                                       newSurvival = if(!is.null(new.survival)){new.survival[indexStrata,,drop=FALSE]}else{new.survival <- matrix(NA)},
                                       firstJumpTime = iid.object$etime1.min[iStrata],
                                       lastSampleTime = iid.object$etime.max[iStrata],
                                       nTau = nTimes,
                                       nNewObs = length(indexStrata),
                                       nSample = object.n,
                                       p = nVar,
                                       exportSE = ("se" %in% export),
                                       exportIF = ("iid" %in% export),
                                       exportIFsum_cumhazard = ("average.iid" %in% export && "cumhazard" %in% type),
                                       exportIFsum_survival = ("average.iid" %in% export && "survival" %in% type),
                                       logTransform = log.transform
                                       )
            if(log.transform){
                if("iid" %in% export){
                    if("cumhazard" %in% type){out$cumhazard.iid[indexStrata,,] <- resCpp$iid}
                    if("survival" %in% type){out$survival.iid[indexStrata,,] <- -resCpp$iid}
                }
                if("se" %in% export){
                    if("cumhazard" %in% type){out$cumhazard.se[indexStrata,] <- resCpp$se}
                    if("survival" %in% type){out$survival.se[indexStrata,] <- resCpp$se}
                }
            }else{
                if("iid" %in% export){
                    if("cumhazard" %in% type){out$cumhazard.iid[indexStrata,,] <-  resCpp$iid}
                    if("survival" %in% type){out$survival.iid[indexStrata,,] <- sliceMultiply_cpp(-resCpp$iid, M = new.survival[indexStrata,,drop=FALSE])}
                }
                if("se" %in% export){
                    if("cumhazard" %in% type){out$cumhazard.se[indexStrata,] <- resCpp$se}
                    if("survival" %in% type){out$survival.se[indexStrata,] <- resCpp$se * new.survival[indexStrata,,drop=FALSE]}
                }                
                if("average.iid" %in% export){ ## average over strata
                    if("cumhazard" %in% type){out$cumhazard.average.iid <- out$cumhazard.average.iid + resCpp$iidsum_cumhazard/n.new}
                    if("survival" %in% type){out$survival.average.iid <- out$survival.average.iid + resCpp$iidsum_survival/n.new}
                }                
            }
        }        

        # }}}
    }else{
        # {{{ other 
        X_IFbeta_mat <- tcrossprod(iid.object$IFbeta, new.LPdata)
        for(iObs in 1:n.new){
            ## print(iObs)
            #NOTE: cannot perfom log transformation if hazard %in% type (error in predictCox)
            iObs.strata <- new.strata[iObs]
            X_IFbeta <- X_IFbeta_mat[,iObs,drop=FALSE]
            if("hazard" %in% type){
                # Evaluate the influence function for the
                # hazard based on the one of the baseline hazard
                if (nVar == 0) {
                    IF_tempo= iid.object$IFhazard[[iObs.strata]]
                }
                else {
                    IF_tempo= (new.eXb[iObs] * (iid.object$IFhazard[[iObs.strata]] + crossprod(t(X_IFbeta),Lambda0$hazard[[iObs.strata]])))
                }
                if("iid" %in% export){
                    out$hazard.iid[iObs,,] <- t(IF_tempo) 
                }    
            }
    
            if("cumhazard" %in% type || "survival" %in% type){
                # Evaluate the influence function for the
                # cumulative hazard based on the one of the cumulative baseline hazard
                if(nVar == 0){
                    IF_tempo <- iid.object$IFcumhazard[[iObs.strata]]
                }else{
                    IF_tempo <- new.eXb[iObs]*(iid.object$IFcumhazard[[iObs.strata]] + crossprod(t(X_IFbeta), Lambda0$cumhazard[[iObs.strata]]))
                }
                if(log.transform){
                    IF_tempo <- rowScale_cpp(IF_tempo, scale = new.cumhazard[iObs,,drop=FALSE])
                    if(any(times<iid.object$etime1.min[iObs.strata])){
                        IF_tempo[,times<iid.object$etime1.min[iObs.strata]] <- 0
                    }
                    #factorTempo <- new.survival[iObs,,drop=FALSE]/(new.survival[iObs,,drop=FALSE]*log(new.survival[iObs,,drop=FALSE]))
                    #IF_tempo.survival <- rowMultiply_cpp(IF_tempo, scale = factorTempo)
                    if("iid" %in% export){
                        if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <-  t(IF_tempo)}
                        if("survival" %in% type){out$survival.iid[iObs,,] <- t(-IF_tempo)}
                    }                    
                    if("se" %in% export){
                        se_tempo <- sqrt(colSums(IF_tempo^2))
                        if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- se_tempo}
                        if("survival" %in% type){out$survival.se[iObs,] <- se_tempo}
                    }
                }else{
                    if("iid" %in% export){
                        if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <-  t(IF_tempo)}
                        if("survival" %in% type){out$survival.iid[iObs,,] <- t(rowMultiply_cpp(-IF_tempo, scale = new.survival[iObs,,drop=FALSE]))}
                    }
                    if("se" %in% export){
                        se_tempo <- sqrt(colSums(IF_tempo^2))
                        if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- se_tempo}
                        if("survival" %in% type){out$survival.se[iObs,] <- se_tempo * new.survival[iObs,,drop=FALSE]}
                    }  
                    if("average.iid" %in% export){ ## average over observations
                        if("cumhazard" %in% type){out$cumhazard.average.iid <- out$cumhazard.average.iid + IF_tempo/n.new}
                        if("survival" %in% type){out$survival.average.iid <- out$survival.average.iid + rowMultiply_cpp(-IF_tempo, scale = new.survival[iObs,,drop=FALSE])/n.new}
                    }
                }
            }
        }

        
        # }}}
    }
      
    ## export
    return(out)
    
}

# }}}

# {{{ Associated functions

# {{{ selectJump

#' @title Evaluate the influence function at selected times
#'
#' @description Evaluate the influence function at selected times
#' @param IF influence function returned by iidCox
#' @param times the times at which the influence function should be assessed
#' @param type can be \code{"hazard"} or/and \code{"cumhazard"}.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @return An object with the same dimensions as IF
#' 
selectJump <- function(IF, times, type){
  
  nStrata <- length(IF$time)
  for(iStrata in 1:nStrata){
      indexJump <- prodlim::sindex(jump.times = IF$time[[iStrata]], eval.times = times)
      
      if(IF$store.iid == "minimal"){          
          IF$calcIFhazard$Elambda0[[iStrata]] <- IF$calcIFhazard$Elambda0[[iStrata]][indexJump,,drop=FALSE]
          IF$calcIFhazard$cumElambda0[[iStrata]] <- IF$calcIFhazard$cumElambda0[[iStrata]][indexJump,,drop=FALSE]
      }else{
          if("hazard" %in% type){
              IFtempo <- matrix(0, nrow = NROW(IF$IFhazard[[iStrata]]), ncol = length(times))
              match.times <- na.omit(match(times, table = IF$time[[iStrata]]))
              if(length(match.times)>0){
                  IFtempo[,times %in% IF$time[[iStrata]]] <- IF$IFhazard[[iStrata]][,match.times,drop=FALSE]
              }
      
              ## name columns
              if(!is.null(colnames(IF$IFhazard[[iStrata]]))){
                  colnames(IFtempo) <- times
              }
      
              ## store
              IF$IFhazard[[iStrata]] <- IFtempo
          }
    
          if("cumhazard" %in% type || "survival" %in% type){
              IF$IFcumhazard[[iStrata]] <- cbind(0,IF$IFcumhazard[[iStrata]])[,indexJump+1,drop = FALSE]
      
              ## name columns
              if(!is.null(colnames(IF$IFcumhazard[[iStrata]]))){
                  colnames(IF$IFcumhazard[[iStrata]]) <- times
              }
          }    
      }
      IF$time[[iStrata]] <- times
  }
  
  return(IF)
  
}

# }}}

#----------------------------------------------------------------------
### calcSeCox.R ends here
