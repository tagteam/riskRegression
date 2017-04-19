# {{{ predictCox
# {{{ header

#' Fast computation of survival probabilities, hazards and cumulative hazards from Cox regression models 
#'
#' Fast routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a \code{survival::coxph} or \code{rms::cph} object
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata A \code{data.frame} or \code{data.table} containing
#'     the values of the predictor variables defining subject specific
#'     predictions. Should have the same structure as the data set
#'     used to fit the \code{object}.
#' @param times Time points at which to evaluate the predictions.
#' @param centered Logical. If \code{TRUE} return prediction at the
#'     mean values of the covariates \code{fit$mean}, if \code{FALSE}
#'     return a prediction for all covariates equal to zero.  in the
#'     linear predictor. Will be ignored if argument \code{newdata} is
#'     used.
#' @param type the type of predicted value. Choices are \itemize{
#'     \item \code{"hazard"} the baseline hazard function when
#'     argument \code{newdata} is not used and the hazard function
#'     when argument \code{newdata} is used.  \item \code{"cumhazard"}
#'     the cumulative baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  \item \code{"survival"}
#'     the survival baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  } Several choices can be
#'     combined in a vector of strings that match (no matter the case)
#'     strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param keep.strata Logical. If \code{TRUE} add the (newdata) strata
#'     to the output. Only if there any.
#' @param keep.times Logical. If \code{TRUE} add the evaluation times
#'     to the output.
#' @param keep.newdata Logical. If \code{TRUE} add the value of the
#'     covariates used to make the prediction in the output list.
#' @param se Logical. If \code{TRUE} add the standard error
#'     corresponding to the output.
#' @param iid Logical. If \code{TRUE} add the influence function
#'     corresponding to the output.
#' @param conf.level Level of confidence.
#' @details Not working with time varying predictor variables or
#'     delayed entry.  The centered argument enables us to reproduce
#'     the results obtained with the \code{basehaz} function from the
#'     survival package.
#'
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' @return 
#' A list with some or all of the following elements:
#' \itemize{
#' \item{times}: the time points at which the other elements are evaluated.
#' \item{hazard}: When argument \code{newdata} is not used the baseline hazard function, otherwise the predicted hazard function. 
#' \item{hazard.se}: The standard errors of the predicted hazard function.
#' \item{cumhazard}: When argument \code{newdata} is not used the cumulative baseline hazard function, otherwise the predicted cumulative hazard function. 
#' \item{cumhazard.se}: The standard errors of the predicted cumulative hazard function.
#' \item{survival}: When argument \code{newdata} is not used the survival probabilities corresponding to the baseline hazard, otherwise the predicted survival probabilities.
#' \item{survival.se}: The standard errors of the predicted survival.
#' \item{strata}: The strata variable.
#' }
#' @examples 
#' library(survival)
#' 
#' set.seed(10)
#' d <- sampleData(40,outcome="survival")
#' nd <- sampleData(4,outcome="survival")
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#' # table(duplicated(d$time))
#'
#' predictCox(fit)
#' predictCox(fit,centered=FALSE,type="hazard")
#' predictCox(fit,centered=TRUE,type="hazard")
#' predictCox(fit, newdata=nd, times=c(3,8),se=TRUE)
#' predictCox(fit, newdata=nd, times = 5,iid=TRUE)
#' 
#' cbind(survival::basehaz(fit),predictCox(fit,type="cumhazard")$cumhazard)
#' 
#' # one strata variable
#' fitS <- coxph(Surv(time,event)~strata(X1)+X2,
#'               data=d, ties="breslow", x = TRUE, y = TRUE)
#' 
#' predictCox(fitS)
#' predictCox(fitS, newdata=nd, times = 1)
#'
#' # two strata variables
#' set.seed(1)
#' d$U=sample(letters[1:5],replace=TRUE,size=NROW(d))
#' d$V=sample(letters[4:10],replace=TRUE,size=NROW(d))
#' nd$U=sample(letters[1:5],replace=TRUE,size=NROW(nd))
#' nd$V=sample(letters[4:10],replace=TRUE,size=NROW(nd))
#' fit2S <- coxph(Surv(time,event)~X1+strata(U)+strata(V)+X2,
#'               data=d, ties="breslow", x = TRUE, y = TRUE)
#'
#' cbind(survival::basehaz(fit2S),predictCox(fit2S,type="cumhazard")$cumhazard)
#' predictCox(fit2S)
#' predictCox(fitS, newdata=nd, times = 3)
#' 
#' 
#' @export

# }}}
predictCox <- function(object,
                       newdata=NULL,
                       times,
                       centered = TRUE,
                       type=c("hazard","cumhazard","survival"),
                       keep.strata = TRUE,
                       keep.times = TRUE,
                       keep.newdata = FALSE,
                       se = FALSE,
                       iid = FALSE,
                       conf.level=0.95){
    status=statusM1=NULL
    
    # {{{ treatment of times and stopping rules

    #### Extract elements from object ####
    # we need:          - the total number of observations, the status and eventtime for each observation
    #                   - the strata corresponding to each observation in the training set
    #                   - the name of each strata
    #                   - the value of the linear predictor for each observation in the training set
    if (se==1L || iid==1L){
        if (missing(newdata)) stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
    }
    infoVar <- CoxVariableName(object)
    is.strata <- infoVar$is.strata
    if (missing(times)) {
        nTimes <- 0
        times <- numeric(0)
    }else{
        nTimes <- length(times)
    }
    needOrder <- (nTimes>0 && is.unsorted(times))
    if (needOrder) {
        oorder.times <- order(order(times))
        times.sorted <- sort(times)
    }else{
        if (nTimes==0)
            times.sorted <- numeric(0)
        else
            times.sorted <- times
    }
    object.n <- CoxN(object)
    object.design <- CoxDesign(object)
    object.status <- object.design[["status"]]
    object.time <- object.design[["stop"]]
    object.strata <- CoxStrata(object, data = NULL, stratavars = infoVar$stratavars)
    object.levelStrata <- levels(object.strata)
    # if we predict the hazard for newdata then there is no need to center the covariates
    object.eXb <- exp(CoxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE})) 
    object.baseEstimator <- CoxBaseEstimator(object) 
    nVar <- length(infoVar$lpvars)
    
    #### checks ####
    if(object.baseEstimator == "exact"){
        stop("Prediction with exact handling of ties is not implemented.\n")
    }
    if(nTimes>0 && any(is.na(times))){
        stop("Missing (NA) values in argument \'times\' are not allowed.\n")
    }
    type <- tolower(type)
    if(!is.null(object$weights)){
        stop("predictCox does not know how to handle Cox models fitted with weights \n")
    }
    if(any(type %in% c("hazard","cumhazard","survival") == FALSE)){
        stop("type can only be \"hazard\", \"cumhazard\" or/and \"survival\" \n") 
    }
    if(any(object.design[,"start"]!=0)){
        stop("do not handle left censoring \n") 
    }
    # if(se == TRUE && ncol(resInfo$modeldata) == 0){
    #   stop("cannot compute the standard error when there are not covariates \n")
    # }
    
    #### Do we want to make prediction for a new dataset ? ####
    # if yes we need to define: - the linear predictor for the new dataset
    #                           - the strata corresponding to each observation in the new dataset
    # (optional for se)         - subset the dataset to the variable involved in the linear predictor and center these variables
    #                           - the influence function             

    # }}}
    # {{{ computation of the baseline hazard
    if(!is.null(newdata)){
        new.n <- NROW(newdata)
        newdata <- as.data.table(newdata)
        new.eXb <- exp(CoxLP(object, data = newdata, center = FALSE))
        
        new.strata <- CoxStrata(object, data = newdata, 
                                sterms = infoVar$sterms, 
                                stratavars = infoVar$stratavars, 
                                levels = object.levelStrata, 
                                stratalevels = infoVar$stratalevels)
        
        new.levelStrata <- levels(new.strata)
    }
    
    #### baseline hazard ####
    nStrata <- length(object.levelStrata)
    if(is.strata){etimes.max <- tapply(object.time, object.strata, max) }else{ etimes.max <- max(object.time) } # last event time
    
    # sort the data
    dt.prepare <- data.table(alltimes = object.time,
                             status = object.status,
                             eXb = object.eXb,
                             strata = as.numeric(object.strata) - 1)
    dt.prepare[, statusM1 := 1-status] # sort by statusM1 such that deaths appear first and then censored events
    data.table::setkeyv(dt.prepare, c("strata", "alltimes", "statusM1"))
    # compute the baseline hazard
    Lambda0 <- baseHaz_cpp(alltimes = dt.prepare$alltimes,
                           status = dt.prepare$status,
                           eXb = dt.prepare$eXb,
                           strata = dt.prepare$strata,
                           nPatients = object.n,
                           nStrata = nStrata,
                           emaxtimes = etimes.max,
                           predtimes = times.sorted,
                           cause = 1,
                           Efron = (object.baseEstimator == "efron"))
    
    # }}}
    #### compute hazard and survival #### 
    if (is.null(newdata)){  
        # {{{ results from the training dataset
        if (!("hazard" %in% type)){ Lambda0$hazard <- NULL } 
        if ("survival" %in% type){  # must be before cumhazard
            Lambda0$survival = exp(-Lambda0$cumhazard)
        }
        if (!("cumhazard" %in% type)){ Lambda0$cumhazard <- NULL } 
        if (keep.times==FALSE){
            Lambda0$time <- NULL
        } 
        if (is.strata == TRUE && keep.strata==1L){ ## rename the strata value with the correct levels
            Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
        }else{
            Lambda0$strata <- NULL
        }
        Lambda0$lastEventTime <- etimes.max
        return(Lambda0)
        # }}}
    } else {
        # {{{ predictions in new dataset
        out <- list()
        if(missing(times) || nTimes==0){
            stop("Time points at which to evaluate the predictions are missing \n")
        }
        
        ## subject specific hazard
        if (is.strata==FALSE){
            if ("hazard" %in% type){
                out$hazard <- (new.eXb %o% Lambda0$hazard)
                if (needOrder) out$hazard <- out$hazard[,oorder.times,drop=0L]
            }
            if ("cumhazard" %in% type || "survival" %in% type){
                cumhazard <- new.eXb %o% Lambda0$cumhazard
                if ("cumhazard" %in% type){
                    if (needOrder)
                        out$cumhazard <- cumhazard[,oorder.times,drop=0L]
                    else
                        out$cumhazard <- cumhazard
                }
                if ("survival" %in% type){
                    out$survival <- exp(-cumhazard)
                    if (needOrder)
                        out$survival <- out$survival[,oorder.times,drop=0L]
                }
            }
            
        }else{ 
            
            ## initialization
            if ("hazard" %in% type){
                out$hazard <- matrix(0, nrow = new.n, ncol = nTimes)
                if(se){
                    out$hazard.se <- matrix(0, nrow = new.n, ncol = nTimes)
                }
            }
            if ("cumhazard" %in% type){
                out$cumhazard <- matrix(NA, nrow = new.n, ncol = nTimes)
                if(se){
                    out$cumhazard.se <- matrix(0, nrow = new.n, ncol = nTimes)
                }
            }
            if ("survival" %in% type){
                out$survival <- matrix(NA, nrow = new.n, ncol = nTimes)
                if(se){
                    out$survival.se <- matrix(0, nrow = new.n, ncol = nTimes)
                }
            }
            if (is.strata == TRUE){ ## rename the strata value with the correct levels
                Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
            }
            
            ## loop across strata
            for(S in new.levelStrata){
                id.S <- Lambda0$strata==S
                newid.S <- new.strata==S
                if ("hazard" %in% type){
                    out$hazard[newid.S,] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S]
                    if (needOrder)
                        out$hazard[newid.S,] <- out$hazard[newid.S,oorder.times,drop=0L]
                }
                if ("cumhazard" %in% type || "survival" %in% type){
                    cumhazard.S <-  new.eXb[newid.S] %o% Lambda0$cumhazard[id.S]
                    if ("cumhazard" %in% type){
                        if (needOrder){
                            out$cumhazard[newid.S,] <- cumhazard.S[,oorder.times,drop=0L]
                        } else{
                            out$cumhazard[newid.S,] <- cumhazard.S
                        }
                    }
                    if ("survival" %in% type){
                        if (needOrder){
                            out$survival[newid.S,] <- exp(-cumhazard.S)[,oorder.times,drop=0L]
                        }else{
                          out$survival[newid.S,] <- exp(-cumhazard.S)
                        }
                    }
                }
            }
        }
        # }}}
        # {{{ standard error
        if(se==1L || iid==1L){ 
            if(nVar > 0){
                # remove response variable
                f.object <- stats::reformulate(attr(stats::terms(CoxFormula(object)),"term.label"),
                                               response = NULL)
                # use prodlim to get the design matrix
                terms.newdata <- stats::terms(f.object, special = CoxSpecialStrata(object), data = newdata)
                new.LPdata <- prodlim::model.design(stats::terms(terms.newdata),
                                                    data = newdata,
                                                    specialsFactor = TRUE,
                                                    dropIntercept = TRUE)$design
                if(NROW(new.LPdata)!=NROW(newdata)){
                    stop("NROW of the design matrix and newdata differ \n",
                         "maybe because newdata contains NA values \n")
                }
            }else{
                new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
            }
            ## influence function 
            iid.object <- object$iid
            
            if(is.null(iid.object)){
                iid.object <- iidCox(object, tauHazard = times.sorted)
            }else{        
                iid.object <- selectJump(iid.object, times = times.sorted, type = type)        
            }
            outSE <- seRobustCox(nTimes = nTimes, type = type,
                                 Lambda0 = Lambda0, iid = iid.object, object.n = object.n, nStrata = nStrata, 
                                 new.eXb = new.eXb, new.LPdata = new.LPdata, new.strata = new.strata, new.survival = out$survival,
                                 nVar = nVar, export = c("iid"[iid==TRUE],"se"[se==TRUE]))
            if(iid == TRUE){
                if ("hazard" %in% type){
                    if (needOrder)
                        out$hazard.iid <- outSE$hazard.iid[,oorder.times,,drop=0L]
                    else
                        out$hazard.iid <- outSE$hazard.iid
                }
                if ("cumhazard" %in% type){
                    if (needOrder)
                        out$cumhazard.iid <- outSE$cumhazard.iid[,oorder.times,,drop=0L]
                    else
                        out$cumhazard.iid <- outSE$cumhazard.iid
                }
                if ("survival" %in% type){
                    if (needOrder)
                        out$survival.iid <- outSE$survival.iid[,oorder.times,,drop=0L]
                    else
                        out$survival.iid <- outSE$survival.iid
                }
            }
            if(se == TRUE){
                zval <- qnorm(1- (1-conf.level)/2, 0,1)
                if ("hazard" %in% type){
                    if (needOrder)
                        out$hazard.se <- outSE$hazard.se[,oorder.times,drop=0L]
                    else
                        out$hazard.se <- outSE$hazard.se
                    
                    out$hazard.lower <- matrix(NA, nrow = NROW(out$hazard.se), ncol = NCOL(out$hazard.se)) # to keep matrix format even when out$hazard contains only one line
                    out$hazard.lower[] <- apply(out$hazard - zval*out$hazard.se,2,pmax,0)
                    out$hazard.upper <- out$hazard + zval*out$hazard.se
                }
                if ("cumhazard" %in% type){
                    if (needOrder)
                        out$cumhazard.se <- outSE$cumhazard.se[,oorder.times,drop=0L]
                    else
                        out$cumhazard.se <- outSE$cumhazard.se
                    
                    out$cumhazard.lower <- matrix(NA, nrow = NROW(out$cumhazard.se), ncol = NCOL(out$cumhazard.se)) # to keep matrix format even when out$cumhazard contains only one line
                    out$cumhazard.lower[] <- apply(out$cumhazard - zval*out$cumhazard.se,2,pmax,0)
                    out$cumhazard.upper <- out$cumhazard + zval*out$cumhazard.se
                }
                if ("survival" %in% type){
                    if (needOrder)
                        out$survival.se <- outSE$survival.se[,oorder.times,drop=0L]
                    else
                        out$survival.se <- outSE$survival.se

                    # to keep matrix format even when out$survival contains only one line
                    out$survival.lower <- out$survival.upper <- matrix(NA, nrow = NROW(out$survival.se), ncol = NCOL(out$survival.se)) 
                    out$survival.lower[] <- apply(out$survival - zval*out$survival.se,2,pmax,0)
                    out$survival.upper[] <- apply(out$survival + zval*out$survival.se,2,pmin,1)
                }
            }
        }
        # }}}
        # {{{ export 
        if (keep.times==TRUE) out <- c(out,list(times=times))
        if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=new.strata))
        out <- c(out,list(lastEventTime=etimes.max,se=se,type=type))
        if( keep.newdata==TRUE){
            out$newdata <- newdata[, CoxCovars(object), with = FALSE]
        }
        class(out) <- "predictCox"
        return(out)
        # }}}
    }
    
}

# }}}

#### Auxiliary functions ####

# {{{ seRobustCox

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
#' @param iid  the value of the influence function returned by \code{iidCox}.
#' @param object.n the number of observations in the dataset used to estimate the object. 
#' @param nStrata the number of strata.
#' @param new.eXb the linear predictor evaluated for the new observations
#' @param new.LPdata the variables involved in the linear predictor for the new observations
#' @param new.strata the strata indicator for the new observations
#' @param new.survival the survival evaluated for the new observations
#' @param nVar the number of variables that form the linear predictor
#' @param export can be "iid" to return the value of the influence function for each observation
#'                      "se" to return the standard error for a given timepoint
#'                      
#'  
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
seRobustCox <- function(nTimes, type, 
                        Lambda0, iid, object.n, nStrata,
                        new.eXb, new.LPdata, new.strata, new.survival,
                        nVar, export){

  n.new <- length(new.eXb)
  
  new.strata <- as.numeric(new.strata)
  
  if(length(Lambda0$strata)==0){
    Lambda0$strata <- rep(1, length(Lambda0$time))
  }else{
    Lambda0$strata <- as.numeric(Lambda0$strata)    
  }
  
  if("hazard" %in% type){Lambda0$hazard <- lapply(1:nStrata,function(s){Lambda0$hazard[Lambda0$strata==s]})}
  if("cumhazard" %in% type || "survival" %in% type){Lambda0$cumhazard <- lapply(1:nStrata,function(s){Lambda0$cumhazard[Lambda0$strata==s]})}
  
  ## main loop
  out <- list()
  if("se" %in% export){
    if("hazard" %in% type){out$hazard.se <- matrix(NA, nrow = n.new, ncol = nTimes)}
    if("cumhazard" %in% type){out$cumhazard.se <- matrix(NA, nrow = n.new, ncol = nTimes)}
    if("survival" %in% type){out$survival.se <- matrix(NA, nrow = n.new, ncol = nTimes)}
  }
  if("iid" %in% export){
    if("hazard" %in% type){out$hazard.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
    if("cumhazard" %in% type){out$cumhazard.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
    if("survival" %in% type){out$survival.iid <- array(NA, dim = c(n.new, nTimes, object.n))}
  }
  
    for(iObs in 1:n.new){
        iObs.strata <- new.strata[iObs]
        X_ICbeta <- iid$ICbeta %*% t(new.LPdata[iObs,,drop=FALSE])
    
      if("hazard" %in% type){      
      IF_tempo <- IClambda2hazard(eXb = new.eXb[iObs],
                                  lambda0 = Lambda0$hazard[[iObs.strata]],
                                  X_ICbeta = X_ICbeta,
                                  IClambda0 = iid$IChazard[[iObs.strata]],
                                  nVar = nVar)
      if("iid" %in% export){
        out$hazard.iid[iObs,,] <- t(IF_tempo)
      }
      if("se" %in% export){
        se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
        out$hazard.se[iObs,] <- se_tempo
      }
    }
    
    if("cumhazard" %in% type || "survival" %in% type){
      IF_tempo <- IClambda2hazard(eXb = new.eXb[iObs],
                                  lambda0 = Lambda0$cumhazard[[iObs.strata]],
                                  X_ICbeta = X_ICbeta,
                                  IClambda0 = iid$ICcumhazard[[iObs.strata]],
                                  nVar = nVar)
      
      if("iid" %in% export){
        if("cumhazard" %in% type){out$cumhazard.iid[iObs,,] <- t(IF_tempo)}
        if("survival" %in% type){out$survival.iid[iObs,,] <- t(rowMultiply_cpp(IF_tempo, scale = new.survival[iObs,,drop=FALSE]))}
      }
      if("se" %in% export){
        se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
        if("cumhazard" %in% type){out$cumhazard.se[iObs,] <- se_tempo}
        if("survival" %in% type){out$survival.se[iObs,] <- se_tempo * new.survival[iObs,,drop=FALSE]}
      }
     
    }
    
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
##' @param X_ICbeta the design matrix times the influence function of beta
##' @param lambda0 the baseline hazard
##' @param IClambda0 the influence function of the baseline hazard 
##' @param nVar the number of variables that form the linear predictor
##' 
##' 
IClambda2hazard <- function(eXb, X_ICbeta, lambda0, IClambda0, nVar){
    if(nVar == 0){
        return(IClambda0)
    }else{
        return(eXb*(IClambda0 + X_ICbeta %*% lambda0))
    }
}

# }}}

# {{{ selectJump

#' @title Evaluate the influence function at selected times
#'
#' @description Evaluate the influence function at selected times
#' @param IC influence function returned by iidCox
#' @param times the times at which the influence function should be assessed
#' @param type can be \code{"hazard"} or/and \code{"cumhazard"}.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @return An object with the same dimensions as IC
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
selectJump <- function(IC, times, type){
  
  nStrata <- length(IC$time)
  for(iStrata in 1:nStrata){
    
    if("hazard" %in% type){
      ICtempo <- matrix(0, nrow = NROW(IC$IChazard[[iStrata]]), ncol = length(times))
      match.times <- na.omit(match(times, table = IC$time[[iStrata]]))
      if(length(match.times)>0){
        ICtempo[,times %in% IC$time[[iStrata]]] <- IC$IChazard[[iStrata]][,match.times,drop=FALSE]
      }
      
      ## name columns
      if(!is.null(colnames(IC$IChazard[[iStrata]]))){
        colnames(ICtempo) <- times
      }
      
      ## store
      IC$IChazard[[iStrata]] <- ICtempo
    }
    
    if("cumhazard" %in% type || "survival" %in% type){
      indexJump <- prodlim::sindex(jump.times = IC$time[[iStrata]], eval.times = times) 
      IC$ICcumhazard[[iStrata]] <- cbind(0,IC$ICcumhazard[[iStrata]])[,indexJump+1,drop = FALSE]
      
      ## name columns
      if(!is.null(colnames(IC$ICcumhazard[[iStrata]]))){
        colnames(IC$ICcumhazard[[iStrata]]) <- times
      }
    }
    
    IC$time[[iStrata]] <- times
  }
  
  
  return(IC)
  
}

# }}}

