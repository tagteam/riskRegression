### calcSeCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (11:46) 
## Version: 
## last-updated: maj 28 2017 (15:16) 
##           By: Brice Ozenne
##     Update #: 128
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
#' Compute the standard error associated to the predictions from Cox regression model using the functional delta method.
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param nTimes the number of time points at which to evaluate the standard errors of the predictions. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param Lambda0 the baseline hazard estimate returned by \code{BaseHazStrata_cpp}.
#' @param object.n the number of observations in the dataset used to estimate the object. 
#' @param nStrata the number of strata.
#' @param new.eXb the linear predictor evaluated for the new observations
#' @param new.LPdata the variables involved in the linear predictor for the new observations
#' @param new.strata the strata indicator for the new observations
#' @param new.cumhazard the cumulative hazard evaluated for the new observations
#' @param new.survival the survival evaluated for the new observations
#' @param nVar the number of variables that form the linear predictor
#' @param logTransform Should the variance/influence function be computed on the log or log(-log) scale
#' @param export can be "iid" to return the value of the influence function for each observation
#'                      "se" to return the standard error for a given timepoint
#'                      
#'  
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
calcSeCox <- function(object, times, nTimes, type, 
                      Lambda0, object.n, object.time, object.eXb, object.strata, nStrata,
                      new.eXb, new.LPdata, new.strata, new.survival, new.cumhazard,
                      nVar, logTransform, export, method.iid){

    # {{{ computation of the influence function
    iid.object <- object$iid            
    if(is.null(iid.object)){
        iid.object <- iidCox(object, tauHazard = times, method.iid = method.iid)
    }else{
        method.iid <- iid.object$method.iid
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
    # }}}

    if(method.iid == "minimal"){
        # {{{ method = minimal
        object.strata <- as.numeric(object.strata)

        if("hazard" %in% type){
            stop("method.iid=\"minimal\" cannot be used to extract the influence function of the hazard \n")
        }
        if(logTransform==FALSE){
            stop("method.iid=\"minimal\" can only be used with logTransform=TRUE \n")
        }
  
        for(iStrata in 1:nStrata){ # iStrata <- 1
            indexStrata <- which(new.strata==iStrata)
 
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
                                       firstJumpTime = iid.object$etime1.min[iStrata],
                                       lastSampleTime = iid.object$etime.max[iStrata],
                                       nTau = nTimes,
                                       nNewObs = length(indexStrata),
                                       nSample = object.n,
                                       p = nVar,
                                       exportSE = ("se" %in% export),
                                       exportIF = ("iid" %in% export)
                                       )

            
            
            if("iid" %in% export){
                if("cumhazard" %in% type){out$cumhazard.iid[indexStrata,,] <- resCpp$iid}
                if("survival" %in% type){out$survival.iid[indexStrata,,] <- resCpp$iid}
            }
            if("se" %in% export){
                if("cumhazard" %in% type){out$cumhazard.se[indexStrata,] <- resCpp$se}
                if("survival" %in% type){out$survival.se[indexStrata,] <- resCpp$se}
            }           

        }        

        # }}}
    }else{
        # {{{ other method
        for(iObs in 1:n.new){
            #NOTE: cannot perfom log transformation if hazard %in% type (error in predictCox)
            
            iObs.strata <- new.strata[iObs]
            X_IFbeta <- iid.object$IFbeta %*% t(new.LPdata[iObs,,drop=FALSE])
        
            if("hazard" %in% type){      
                IF_tempo <- IFlambda2hazard(eXb = new.eXb[iObs],
                                            lambda0 = Lambda0$hazard[[iObs.strata]],
                                            X_IFbeta = X_IFbeta,
                                            IFlambda0 = iid.object$IFhazard[[iObs.strata]],
                                            nVar = nVar)
                if("iid" %in% export){
                    out$hazard.iid[iObs,,] <- t(IF_tempo) 
                }    
            }
    
            if("cumhazard" %in% type || "survival" %in% type){
                IF_tempo <- IFlambda2hazard(eXb = new.eXb[iObs],
                                            lambda0 = Lambda0$cumhazard[[iObs.strata]],
                                            X_IFbeta = X_IFbeta,
                                            IFlambda0 = iid.object$IFcumhazard[[iObs.strata]],
                                            nVar = nVar)
                
                if(logTransform){
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
                        se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
                        if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- se_tempo}
                        if("survival" %in% type){out$survival.se[iObs,] <- se_tempo}
                    }
                }else{
                    if("iid" %in% export){
                        if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <-  t(IF_tempo)}
                        if("survival" %in% type){out$survival.iid[iObs,,] <- t(rowMultiply_cpp(IF_tempo, scale = new.survival[iObs,,drop=FALSE]))}
                    }
                    if("se" %in% export){
                        se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
                        if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- se_tempo}
                        if("survival" %in% type){out$survival.se[iObs,] <- se_tempo * new.survival[iObs,,drop=FALSE]}
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

# {{{ lambda2hazard

##' @title Evaluate the influence function for the hazard based on the one of the baseline hazard##' 
##' @description Evaluate the influence function for the hazard based on the one of the baseline hazard
##'
##' @param eXb the linear predictor
##' @param X_IFbeta the design matrix times the influence function of beta
##' @param lambda0 the baseline hazard
##' @param IFlambda0 the influence function of the baseline hazard 
##' @param nVar the number of variables that form the linear predictor
##' 
##' 
IFlambda2hazard <- function(eXb, X_IFbeta, lambda0, IFlambda0, nVar){
    if(nVar == 0){
        return(IFlambda0)
    }else{
        return(eXb*(IFlambda0 + X_IFbeta %*% lambda0))
    }
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
#' @examples 
#' \dontrun{
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow", x = TRUE, y = TRUE)
#' 
#' IFall <- iidCox(fit)
#' selectJump(IFall, times = 1:2, type = "cumhazard") 
#'  
#' }
selectJump <- function(IF, times, type){
  
  nStrata <- length(IF$time)
  for(iStrata in 1:nStrata){
      indexJump <- prodlim::sindex(jump.times = IF$time[[iStrata]], eval.times = times)
      
      if(IF$method.iid == "minimal"){          
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

# }}}
#----------------------------------------------------------------------
### calcSeCox.R ends here
