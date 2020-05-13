### calcSeCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (11:46) 
## Version: 
## last-updated: maj 13 2020 (09:39) 
##           By: Brice Ozenne
##     Update #: 707
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


## * calcSeCox (documentation)
#' @title Computation of standard errors for predictions
#' @description Compute the standard error associated to the predictions from Cox regression model
#' using a first order von Mises expansion of the functional (cumulative hazard or survival).
#' @name calcSeCox
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param times Vector of times at which to return the estimated
#'      hazard/survival.
#' @param nTimes the length of the argument \code{times}. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param diag [logical] when \code{FALSE} the hazard/cumlative hazard/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param Lambda0 the baseline hazard estimate returned by \code{BaseHazStrata_cpp}.
#' @param object.n the number of observations in the dataset used to estimate the object. 
#' @param object.time the time to event of the observations used to estimate the object.
#' @param object.eXb the exponential of the linear predictor relative to the observations used to estimate the object. 
#' @param object.strata the strata index of the observations used to estimate the object.
#' @param nStrata the number of strata.
#' @param new.n the number of observations for which the prediction was performed.
#' @param new.eXb the linear predictor evaluated for the new observations.
#' @param new.LPdata the variables involved in the linear predictor for the new observations.
#' @param new.strata the strata indicator for the new observations.
#' @param new.survival the survival evaluated for the new observations.
#' @param nVar the number of variables that form the linear predictor.
#' @param export can be "iid" to return the value of the influence function for each observation.
#'                      "se" to return the standard error for a given timepoint.
#'                      
#' @details Can also return the estimated influence function for the cumulative hazard function and survival probabilities
#'  the sum over the observations of the estimated influence function.
#'
#' \code{store.iid="full"} compute the influence function for each observation at each time in the argument \code{times}
#' before computing the standard error / influence functions.
#' \code{store.iid="minimal"} recompute for each subject specific prediction the influence function for the baseline hazard.
#' This avoid to store all the influence functions but may lead to repeated evaluation of the influence function.
#' This solution is therefore more efficient in memory usage but may not be in terms of computation time.
#' 
# #' @inheritParams  predict.CauseSpecificCox
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
#'

## * calcSeCox (code)
#' @rdname calcSeCox
calcSeCox <- function(object, times, nTimes, type, diag,
                      Lambda0, object.n, object.time, object.eXb, object.strata, nStrata,
                      new.n, new.eXb, new.LPdata, new.strata, new.survival, 
                      nVar, export){

    ## ** Prepare arguments
    if(diag){nTimes <- 1}
    new.strata <- as.numeric(new.strata)
    
    out <- list()
    if("se" %in% export){
        if("cumhazard" %in% type){out$cumhazard.se <- matrix(NA, nrow = new.n, ncol = nTimes)}
        if("survival" %in% type){out$survival.se <- matrix(NA, nrow = new.n, ncol = nTimes)}
    }
    if("iid" %in% export){
        if("hazard" %in% type){out$hazard.iid <- array(NA, dim = c(new.n, nTimes, object.n))}
        if("cumhazard" %in% type){out$cumhazard.iid <- array(NA, dim = c(new.n, nTimes, object.n))}
        if("survival" %in% type){out$survival.iid <- array(NA, dim = c(new.n, nTimes, object.n))}
    }
    if("average.iid" %in% export){
        if("hazard" %in% type){out$hazard.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)}
        if("cumhazard" %in% type){out$cumhazard.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)}
        if("survival" %in% type){out$survival.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)}
    }

    if(nVar>0){
        X_IFbeta_mat <- tcrossprod(object$iid$IFbeta, new.LPdata)
    }
        
    if(diag){

         for(iStrata in 1:nStrata){ ## iStrata <- 1
            indexStrata <- which(new.strata==iStrata)                
            if(length(indexStrata)==0){next}
            iPrevalence <- length(indexStrata)/new.n

            ## compute iid
            if("hazard" %in% type){
                if(nVar==0){
                    iIFhazard <- object$iid$IFhazard[[iStrata]][,indexStrata,drop=FALSE]
                }else{
                    iIFhazard <- rowMultiply_cpp(object$iid$IFhazard[[iStrata]][,indexStrata,drop=FALSE] + rowMultiply_cpp(X_IFbeta_mat[,indexStrata], scale = Lambda0$hazard[[iStrata]][indexStrata]),
                                                 scale = new.eXb[indexStrata])
                }
            }

            if("cumhazard" %in% type || "survival" %in% type){
                if(nVar==0){
                    iIFcumhazard <- object$iid$IFcumhazard[[iStrata]][,indexStrata,drop=FALSE]
                }else{
                    iIFcumhazard <- rowMultiply_cpp(object$iid$IFcumhazard[[iStrata]][,indexStrata,drop=FALSE] + rowMultiply_cpp(X_IFbeta_mat[,indexStrata,drop=FALSE], scale = Lambda0$cumhazard[[iStrata]][indexStrata]),
                                                    scale = new.eXb[indexStrata])
                }
                if("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)){
                    iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[indexStrata,])
                }
            }

            ## export
            if("iid" %in% export){
                if("hazard" %in% type){out$hazard.iid[indexStrata,1,] <- t(iIFhazard)}
                if("cumhazard" %in% type){out$cumhazard.iid[indexStrata,1,] <- t(iIFcumhazard)}
                if("survival" %in% type){out$survival.iid[indexStrata,1,] <- t(iIFsurvival)}
            }
            if("se" %in% export){
                iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
                if("cumhazard" %in% type){out$cumhazard.se[indexStrata,1] <- iSEcumhazard}
                if("survival" %in% type){out$survival.se[indexStrata,1] <- iSEcumhazard * new.survival[indexStrata,1]}
            }
            if("average.iid" %in% export){
                if("hazard" %in% type){out$hazard.average.iid[,1] <- out$hazard.average.iid[,1] + rowSums(iIFhazard)/new.n}
                if("cumhazard" %in% type){out$cumhazard.average.iid[,1] <- out$cumhazard.average.iid[,1] + rowSums(iIFcumhazard)/new.n}
                if("survival" %in% type){out$survival.average.iid[,1] <- out$survival.average.iid[,1] + rowSums(iIFsurvival)/new.n}
            }
         }
        
    }else if(nVar==0){

        for(iStrata in 1:nStrata){ ## iStrata <- 1
            indexStrata <- which(new.strata==iStrata)                
            if(length(indexStrata)==0){next}
            iPrevalence <- length(indexStrata)/new.n

            ## compute iid
            if("hazard" %in% type){
                iIFhazard <- object$iid$IFhazard[[iStrata]]
                tiIFhazard <- t(iIFhazard)
            }
            
            if("cumhazard" %in% type || "survival" %in% type){
                iIFcumhazard <- object$iid$IFcumhazard[[iStrata]]
                tiIFcumhazard <- t(iIFcumhazard)

                if("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)){
                    iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[indexStrata[1],]) ## everybody has the same survival
                    tiIFsurvival <- t(iIFsurvival)
                }
            }
                
            ## export
            if("se" %in% export){
                iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
                if("survival" %in% type){
                    iSEsurvival <- iSEcumhazard * new.survival[indexStrata[1],] ## everybody has the same survival
                }
            }

            for(iObs in indexStrata){ ## iObs <- 1
                if("iid" %in% export){
                    if("hazard" %in% type){out$hazard.iid[iObs,,] <- tiIFhazard}
                    if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <- tiIFcumhazard}
                    if("survival" %in% type){out$survival.iid[iObs,,] <- tiIFsurvival}
                }
                if("se" %in% export){
                    if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- iSEcumhazard}
                    if("survival" %in% type){out$survival.se[iObs,] <- iSEsurvival}
                }                        
            }
            if("average.iid" %in% export){
                if("hazard" %in% type){out$hazard.average.iid <- out$hazard.average.iid + iIFhazard * iPrevalence}
                if("cumhazard" %in% type){out$cumhazard.average.iid <- out$cumhazard.average.iid + iIFcumhazard * iPrevalence}
                if("survival" %in% type){out$survival.average.iid <- out$survival.average.iid + iIFsurvival * iPrevalence}
            }
        }

        }else{ ## !diag && nVar>0

            for(iObs in 1:new.n){ ## iObs <- 1
                iObs.strata <- new.strata[iObs]

                ## compute iid
                if("hazard" %in% type){
                    iIFhazard <- (new.eXb[iObs] * (object$iid$IFhazard[[iObs.strata]] + crossprod(t(X_IFbeta_mat[,iObs,drop=FALSE]),Lambda0$hazard[[iObs.strata]])))
                }             
                if("cumhazard" %in% type || "survival" %in% type){
                    iIFcumhazard <- new.eXb[iObs] * (object$iid$IFcumhazard[[iObs.strata]] + crossprod(t(X_IFbeta_mat[,iObs,drop=FALSE]), Lambda0$cumhazard[[iObs.strata]]))
                }
                if("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)){
                    iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[iObs,])
                }
                
                ## export                    
                if("iid" %in% export){
                    if("hazard" %in% type){out$hazard.iid[iObs,,] <- t(iIFhazard)}
                    if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <- t(iIFcumhazard)}
                    if("survival" %in% type){out$survival.iid[iObs,,] <- t(iIFsurvival)}
                }
                if("se" %in% export){
                    iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
                    if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- iSEcumhazard}
                    if("survival" %in% type){out$survival.se[iObs,] <- iSEcumhazard * new.survival[iObs,,drop=FALSE]}
                }
                if("average.iid" %in% export){  ## average over observations
                    if("hazard" %in% type){out$hazard.average.iid <- out$hazard.average.iid + iIFhazard/new.n}
                    if("cumhazard" %in% type){out$cumhazard.average.iid <- out$cumhazard.average.iid + iIFcumhazard/new.n}
                    if("survival" %in% type){out$survival.average.iid <- out$survival.average.iid + iIFsurvival/new.n}
                }
            }
        }

    ## export
    return(out)
}


## * calcAiidCox (code)
## Average iid 
#' @rdname calcSeCox
calcAiidCox <- function(object, times, nTimes, type, diag,
                        Lambda0, object.n, object.time, object.eXb, object.strata, nStrata,
                        new.n, new.eXb, new.LPdata, new.strata, new.survival, 
                        nVar, factor){


    if(diag){nTimes <- 1}

    ## ** prepare strata 
    new.Ustrata <- sort(unique(new.strata))
    new.nStrata <- length(new.Ustrata)
            
    new.indexStrata <- lapply(new.Ustrata, function(iStrata){
        which(new.strata==iStrata) - 1
    })
    new.prevStrata <- sapply(new.indexStrata, length)/new.n           

    ## ** normalize arguments for C++
    attr(new.LPdata,"levels") <- NULL
    if(is.null(new.survival)){
        new.survival <- matrix()
    }

    if(is.null(factor)){
        rm.list <- TRUE
        factor <- list(matrix(1, nrow = new.n, ncol = nTimes))
    }else{
        rm.list <- FALSE                
    }

    ## ** C++
    if("hazard" %in% type){
        outRcpp.hazard <- calcAIFsurv_cpp(ls_IFcumhazard = object$iid$IFhazard[new.Ustrata], 
                                          IFbeta = object$iid$IFbeta,
                                          cumhazard0 = Lambda0$hazard[new.Ustrata],
                                          survival = matrix(0),
                                          eXb = new.eXb,
                                          X = new.LPdata,
                                          prevStrata = new.prevStrata,
                                          ls_indexStrata = new.indexStrata,
                                          factor = factor,
                                          nTimes = nTimes,
                                          nObs = object.n,
                                          nStrata = new.nStrata,
                                          nVar = nVar,
                                          diag = diag,
                                          exportCumHazard = TRUE,
                                          exportSurvival = FALSE)
    }
    if(("cumhazard" %in% type) || ("survival" %in% type)){
        outRcpp.cumhazard <- calcAIFsurv_cpp(ls_IFcumhazard = object$iid$IFcumhazard[new.Ustrata], 
                                             IFbeta = object$iid$IFbeta,
                                             cumhazard0 = Lambda0$cumhazard[new.Ustrata],
                                             survival = new.survival,
                                             eXb = new.eXb,
                                             X = new.LPdata,
                                             prevStrata = new.prevStrata,
                                             ls_indexStrata = new.indexStrata,
                                             factor = factor,
                                             nTimes = nTimes,
                                             nObs = object.n,
                                             nStrata = new.nStrata,
                                             nVar = nVar,
                                             diag = diag,
                                             exportCumHazard = ("cumhazard" %in% type),
                                             exportSurvival = ("survival" %in% type))
    }
        
    ## ** reshape
    out <- list()
    if("hazard" %in% type){
        if(rm.list){
            out$hazard.average.iid <- matrix(outRcpp.hazard[[1]][[1]], nrow = object.n, ncol = nTimes)
        }else{
            out$hazard.average.iid <- lapply(outRcpp.hazard[[1]], function(iMat){matrix(iMat, nrow = object.n, ncol = nTimes)})
        }
    }
    if("cumhazard" %in% type){
        if(rm.list){
            out$cumhazard.average.iid <- matrix(outRcpp.cumhazard[[1]][[1]], nrow = object.n, ncol = nTimes)
        }else{
            out$cumhazard.average.iid <- lapply(outRcpp.cumhazard[[1]], function(iMat){matrix(iMat, nrow = object.n, ncol = nTimes)})
        }
    }
    if("survival" %in% type){
        if(rm.list){
            out$survival.average.iid <- matrix(outRcpp.cumhazard[[2]][[1]], nrow = object.n, ncol = nTimes)
        }else{
            out$survival.average.iid <- lapply(outRcpp.cumhazard[[2]], function(iMat){matrix(iMat, nrow = object.n, ncol = nTimes)})
        }
    }
    
    return(out)
}

## * selectJump
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
      
      if(IF$store.iid == "minimal"){          
          indexJump <- prodlim::sindex(jump.times = IF$time[[iStrata]], eval.times = times)
          IF$calcIFhazard$Elambda0[[iStrata]] <- IF$calcIFhazard$Elambda0[[iStrata]][indexJump,,drop=FALSE]
          IF$calcIFhazard$cumElambda0[[iStrata]] <- IF$calcIFhazard$cumElambda0[[iStrata]][indexJump,,drop=FALSE]
      }else{
          if("hazard" %in% type){
              match.times <- match(times, table = IF$time[[iStrata]])
              match.times[is.na(match.times)] <- 0
              if(any(times > IF$etime.max[[iStrata]])){
                  match.times[times > IF$etime.max[[iStrata]]] <- NA
              }
              IF$IFhazard[[iStrata]] <- subsetIndex(IF$IFhazard[[iStrata]], index = match.times, default = 0, col = TRUE)
          }
          if("cumhazard" %in% type || "survival" %in% type){
              indexJump <- prodlim::sindex(jump.times = IF$time[[iStrata]], eval.times = times)
              if(any(times > IF$etime.max[[iStrata]])){
                  indexJump[times > IF$etime.max[[iStrata]]] <- NA
              }
              IF$IFcumhazard[[iStrata]] <- subsetIndex(IF$IFcumhazard[[iStrata]], index = indexJump, default = 0, col = TRUE)
          }    
      }
      IF$time[[iStrata]] <- times
  }
  
  return(IF)
  
}

#----------------------------------------------------------------------
### calcSeCox.R ends here
