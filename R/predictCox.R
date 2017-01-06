#' Fast prediction of survival, hazard and cumulative hazard from Cox regression model 
#'
#' Fast routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a coxph or cph object
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata  A data frame or table containing the values of the
#'     predictor variables defining subject specific predictions. Should have
#'     the same structure as the data set used to fit the \code{object}.
#' @param times Time points at which to evaluate the predictions. 
#' @param centered If TRUE remove the centering factor used by \code{coxph}
#'     in the linear predictor. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}
#' @param keep.strata Logical. If \code{TRUE} add the (newdata) strata to the output. Only if there any. 
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output. 
#' @param keep.lastEventTime Logical. If \code{TRUE} add the time at which the last event occured in each strata to the output list. 
#' @param format how to export the baseline hazard. Can be \code{"data.frame"}, \code{"data.table"} or \code{"list"}.
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. 
#' @param iid the value of the influence function obtained with the function \code{iidCox}.
#' @details Not working with time varying predictor variables or
#'     delayed entry.
#' The centered argument enables to reproduce results obtained with the \code{basehaz} function from the survival package.
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' @return A list optionally containing the time, the strata (if any), the hazard, the
#'         cumulative hazard and survival probabilities.
#' @examples 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' nd <- SimSurv(10)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow", x = TRUE, y = TRUE)
#' # table(duplicated(d$time))
#' 
#' predictCox(fit)
#' predictCox(fit, newdata=nd, times = 5)
#' cbind(survival::basehaz(fit),predictCox(fit,type="cumHazard"))
#' 
#' # one strata variable
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow", x = TRUE, y = TRUE)
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
#' fit2S <- coxph(Surv(time,status)~X1+strata(U)+strata(V)+X2,data=d, ties="breslow", x = TRUE, y = TRUE)
#'
#' cbind(survival::basehaz(fit2S),predictCox(fit2S,type="cumHazard"))
#' predictCox(fit2S)
#' predictCox(fitS, newdata=nd, times = 3)
#' 
#' 
#' @export
predictCox <- function(object,
                       newdata=NULL,
                       times,
                       centered = TRUE,
                       type=c("hazard","cumHazard","survival"),
                       keep.strata = TRUE,
                       keep.times = TRUE,
                       keep.lastEventTime = FALSE,
                       se = FALSE,
                       iid = NULL,
                       format = "data.frame"){ 
  
  #### extract elements from object ####
  # we need:          - the total number of observations, the status and eventtime for each observation
  #                   - the strata corresponding to each observation in the training set
  #                   - the name of each strata
  #                   - the value of the linear predictor for each observation in the training set
  infoVar <- CoxVariableName(object)
  is.strata <- infoVar$is.strata
  
  
  object.n <- CoxN(object)
  object.design <- CoxDesign(object)
  object.status <- object.design[["status"]]
  object.time <- object.design[["stop"]]
  object.strata <- CoxStrata(object, stratavars = infoVar$stratavars)
  object.levelStrata <- levels(object.strata)
  object.eXb <- exp(CoxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE})) # if we use the linear predictor then no need to center since the centering term will cancel between the linear predictor and the baseline hazard
  object.baseEstimator <- CoxBaseEstimator(object) 
  nVar <- length(infoVar$lpvars)
  
  #### checks ####
  if(object.baseEstimator == "exact"){
    stop("Prediction with exact correction for ties is not implemented \n")
  }
  if(!missing(times) && any(is.na(times))){
    stop("NA values in argument \'times\' \n")
  }
  if(!is.null(newdata) || "eXb" %in% type || "newstrata" %in% type){ # if the baseline hazard is exported
    if(format %in% c("data.frame","data.table","list") == FALSE){
      stop("format can only be \'data.frame\', \'data.table\' , or \'list\' \n")
    }
    if(format %in% c("data.frame","data.table") && keep.lastEventTime){
      stop("format must be \'list\' when \'keep.lastEventTime\' equals TRUE \n")
    }
  }
  if(!is.null(object$weights)){
    stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if(any(type %in% c("hazard","cumHazard","survival") == FALSE)){
    stop("type can only be \"hazard\", \"cumHazard\" or/and \"survival\" \n") 
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
 
  Lambda0 <- baseHaz_cpp(alltimes = object.time,
                         status = object.status,
                         eXb = object.eXb,
                         strata = as.numeric(object.strata) - 1,
                         se = FALSE, # to be ignored
                         data = matrix(0, ncol = 1, nrow = object.n), # to be ignored
                         nVar = 1, # to be ignored
                         nPatients = object.n,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = if(missing(times)){numeric(0)}else{sort(times)},
                         cause = 1,
                         Efron = (object.baseEstimator == "efron"))
  
  Lambda0$Xbar <- NULL
  Lambda0$XbarCumSumRes <- NULL
  Lambda0$se.hazard <- NULL
  Lambda0$se.cumHazard <- NULL
  
  if (is.strata == TRUE){ ## rename the strata value with the correct levels
    Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
  }
  
  #### compute hazard and survival #### 
  if (is.null(newdata)){  ## on the training dataset
    
    if ("hazard" %in% type == FALSE){ Lambda0$hazard <- NULL } 
    
    if ("survival" %in% type){  # must be before cumHazard
      Lambda0$survival = exp(-Lambda0$cumHazard)
    } 
    
    if ("cumHazard" %in% type == FALSE){ Lambda0$cumHazard <- NULL } 
    
    
    if (keep.times==FALSE){
      Lambda0$time <- NULL
    } 
    
    if(is.strata == FALSE || keep.strata==FALSE){
      Lambda0$strata <- NULL
    }
    
    if( keep.lastEventTime==TRUE){
      Lambda0$lastEventTime <- etimes.max
    } 
    
    return(do.call(paste0("as.",format), list(Lambda0)))
    
  } else { ## on a new dataset
    out <- list()
    
    if(missing(times) || length(times)==0){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    nTimes <- length(times)
    
    ## subject specific hazard
    if (is.strata==FALSE){
      
      if ("hazard" %in% type){out$hazard <- (new.eXb %o% Lambda0$hazard)}
      if ("cumHazard" %in% type || "survival" %in% type){
        cumHazard <- new.eXb %o% Lambda0$cumHazard
        if ("cumHazard" %in% type){out$cumHazard <- cumHazard}
        if ("survival" %in% type){out$survival <- exp(-cumHazard)}
      }
      
    }else{ 
      
      ## initialization
      if ("hazard" %in% type){
        out$hazard <- matrix(0, nrow = new.n, ncol = nTimes)
        if(se){out$hazard.se <- matrix(0, nrow = new.n, ncol = nTimes)}
      }
      if ("cumHazard" %in% type){
        out$cumHazard <- matrix(NA, nrow = new.n, ncol = nTimes)
        if(se){out$cumHazard.se <- matrix(0, nrow = new.n, ncol = nTimes)}
      }
      if ("survival" %in% type){
        out$survival <- matrix(NA, nrow = new.n, ncol = nTimes)
        if(se){out$survival.se <- matrix(0, nrow = new.n, ncol = nTimes)}
      }
      
      ## loop across strata
      for(S in new.levelStrata){
        id.S <- Lambda0$strata==S
        newid.S <- new.strata==S
        
        if ("hazard" %in% type){out$hazard[newid.S,] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S]}
        if ("cumHazard" %in% type || "survival" %in% type){
          cumHazard.S <-  new.eXb[newid.S] %o% Lambda0$cumHazard[id.S]
          if ("cumHazard" %in% type){out$cumHazard[newid.S,] <- cumHazard.S}
          if ("survival" %in% type){out$survival[newid.S,] <- exp(-cumHazard.S)}
        }
      }
    }
    
    #### standard error ####
    if(se){ 
      
      if(nVar > 0){
        new.LPdata <- model.matrix(object, newdata)
      }else{
        new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
      }
      
      ## influence function 
      if(is.null(iid)){
        if("hazard" %in% type){
          iid <- iidCox(object)
          iid$IChazard <- calcIChazard(iid$ICcumHazard)
          iid <- selectJump(iid, times = times, type = type)
        }else{
          iid <- iidCox(object, tauHazard = times)
        }
        
      }else{
        
        if("hazard" %in% type){
          iid$IChazard <- calcIChazard(iid$ICcumHazard)
        }
        iid <- selectJump(iid, times = times, type = type)
        
      }
      
      outSE <- seRobustCox(object,  times = times, type = type,
                           Lambda0 = Lambda0, iid = iid, nStrata = nStrata,
                           new.eXb = new.eXb, new.LPdata = new.LPdata, new.strata = new.strata, new.survival = out$survival)
      
      if ("hazard" %in% type){out$hazard.se <- outSE$hazard.se}
      if ("cumHazard" %in% type){out$cumHazard.se <- outSE$cumHazard.se}
      if ("survival" %in% type){out$survival.se <- outSE$survival.se}
    }
    
    #### export ####
    ## if necessary reorder columns according to time
    if(any(order(times) != 1:length(times))){
      if ("hazard" %in% type){
        out$hazard <- out$hazard[,order(order(times)), drop = FALSE]
        if(se){out$hazard.se <- outSE$hazard.se[,order(order(times)), drop = FALSE]}
      }
      if ("cumHazard" %in% type){
        out$cumHazard <- out$cumHazard[,order(order(times)), drop = FALSE]
        if(se){out$cumHazard.se <- outSE$cumHazard.se[,order(order(times)), drop = FALSE]}
      }
      if ("survival" %in% type){
        out$survival <- out$survival[,order(order(times)), drop = FALSE]
        if(se){out$survival.se <- outSE$survival.se[,order(order(times)), drop = FALSE]}
      }
    }
    if (keep.times==TRUE) out <- c(out,list(times=times))
    if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=new.strata))
    if( keep.lastEventTime==TRUE) out <- c(out,list(lastEventTime=etimes.max))
    
    return(out)
  }
  
  
  
}

#' Computation of standard errors for predictions
#'
#' Compute the standard error associated to the predictions from Cox regression model using the functional delta method.
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param times Time points at which to evaluate the standard errors of the predictions. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param Lambda0 the baseline hazard estimate returned by \code{BaseHazStrata_cpp}.
#' @param iid  the value of the influence function returned by \code{iidCox}.
#' @param nStrata	the number of strata.
#' @param new.eXb the linear predictor evaluated for the new observations
#' @param new.LPdata the variables involved in the linear predictor for the new observations
#' @param new.strata the strata indicator for the new observations
#' @param new.survival the survival evaluated for the new observations
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
seRobustCox <- function(object, times, type, 
                        Lambda0, iid, nStrata,
                        new.eXb, new.LPdata, new.strata, new.survival){
  
  n.times <- length(times)
  n.new <- length(new.eXb)
  
  new.strata <- as.numeric(new.strata)
  
  if(length(Lambda0$strata)==0){
    Lambda0$strata <- rep(1, length(Lambda0[[type[1]]]))
  }else{
    Lambda0$strata <- as.numeric(Lambda0$strata)    
  }
  if("hazard" %in% type){Lambda0$hazard <- lapply(1:nStrata,function(s){Lambda0$hazard[Lambda0$strata==s]})}
  if("cumHazard" %in% type){Lambda0$cumHazard <- lapply(1:nStrata,function(s){Lambda0$cumHazard[Lambda0$strata==s]})}
  
  ## main loop
  out <- list()
  if("hazard" %in% type){
    out$hazard.se <- matrix(NA, nrow = n.new, ncol = n.times)
  }
  if("cumHazard" %in% type){
    out$cumHazard.se <- matrix(NA, nrow = n.new, ncol = n.times)
  }
  if("survival" %in% type){
    out$survival.se <- matrix(NA, nrow = n.new, ncol = n.times)
  }
  
  for(iObs in 1:n.new){
    iObs.strata <- new.strata[iObs]
    X_ICbeta <- iid$ICbeta %*% t(new.LPdata[iObs,,drop=FALSE])
    
    # browser()
    if("hazard" %in% type){
      IF_tempo <- new.eXb[iObs]*(iid$IChazard[[iObs.strata]] + X_ICbeta %*% Lambda0$hazard[[iObs.strata]])
      se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
      
      out$hazard.se[iObs,] <- se_tempo
    }
    
    if("cumHazard" %in% type || "survival" %in% type){
      IF_tempo <- new.eXb[iObs]*(iid$ICcumHazard[[iObs.strata]] + X_ICbeta %*% Lambda0$cumHazard[[iObs.strata]])
      se_tempo <- sqrt(apply(IF_tempo^2,2,sum))
      
      if("cumHazard" %in% type){
        out$cumHazard.se[iObs,] <- se_tempo
      }
      
      if("survival" %in% type){
        out$survival.se[iObs,] <- se_tempo * new.survival[iObs,,drop=FALSE]
      }
    }
    
  }
  
  ## export
  return(out)
}


#' @title Evaluate the influence function for the hazard functional
#'
#' @description  Evaluate the influence function for the hazard functional using the one of the cumulative hazard
#' @param ICcumHazard influence function of the cumulative hazard
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @return An object with the same dimensions as ICcumHazard, i.e. a list with one element per strata each element being a matrix observation*time.
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
#' IFhaz <-  calcIChazard(IFall$ICLambda0)
#' 
#' cbind(IFhaz[[1]][,3], IFall$ICLambda0[[1]][,3]-IFall$ICLambda0[[1]][,2])
#'  
#' }
calcIChazard <- function(ICcumHazard){
  
  nStrata <- length(ICcumHazard)
  IChazard <- lapply(1:nStrata, function(x){NULL})
  
  for(iStrata in 1:nStrata){
    IChazard[[iStrata]] <- ICcumHazard[[iStrata]][,1]
    if(ncol(ICcumHazard[[iStrata]])>1){
      IChazard[[iStrata]] <- cbind(IChazard[[iStrata]], t(apply(ICcumHazard[[iStrata]], 1, diff)))
    }
  } 
  colnames(IChazard) <- colnames(ICcumHazard) 
  
  return(IChazard)
  
}

#' @title Evaluate the influence function at selected times
#'
#' @description Evaluate the influence function at selected times
#' @param IC influence function returned by iidCox
#' @param times the times at which the influence function should be assessed
#' @param type can be \code{"hazard"} or/and \code{"cumHazard"}.
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
#' selectJump(IFall, times = 1:2, type = "cumHazard") 
#'  
#' }
selectJump <- function(IC, times, type){
  
  nStrata <- length(IC$time)
  for(iStrata in 1:nStrata){
    indexJump <- prodlim::sindex(jump.times = IC$time[[iStrata]], eval.times = times) 
    
    if("hazard" %in% type){
      IC$IChazard[[iStrata]] <- cbind(0,IC$IChazard[[iStrata]])[,indexJump+1,drop = FALSE]
    }
    if("cumHazard" %in% type || "survival" %in% type){
      IC$ICcumHazard[[iStrata]] <- cbind(0,IC$ICcumHazard[[iStrata]])[,indexJump+1,drop = FALSE]
    }
    IC$time[[iStrata]] <- times
  }
  
  return(IC)
  
}

#' #' Computation of standard errors for predictions
#' #'
#' #' Compute the standard error associated to the predictions from Cox regression model using the functional delta method.
#' #' @param object The fitted Cox regression model object either
#' #'     obtained with \code{coxph} (survival package) or \code{cph}
#' #'     (rms package).
#' #' @param newdata  A data frame or table containing the values of the
#' #'     predictor variables defining subject specific predictions. Should have
#' #'     the same structure as the data set used to fit the \code{object}.
#' #' @param times Time points at which to evaluate the predictions. 
#' #' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' #' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' #' @param Lambda0 the baseline hazard estimate returned by BaseHazStrata_cpp.
#' #' @param eXb exponential (covariates * coefficients)
#' #' @param survival the predicted survival (only necessary when type contains survival)
#' #' @param subset.Lambda0 an index giving the subset of observations to be used in Lambda0
#' #' 
#' #' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' #' 
#' #' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
#' seCox <- function(object, newdata, times, type, Lambda0, eXb, survival, subset.Lambda0 = 1:length(Lambda0$time)){
#'   
#'   ## WARNING newdata must be centered
#'   n.times <- length(times)
#'   n.newdata <- NROW(newdata)
#'   etimes <- Lambda0$time[subset.Lambda0]
#'   newdata <- as.matrix(newdata)
#'   
#'   ##
#'   if("hazard" %in% type){
#'     se.hazard0.tindex <- Lambda0$se.hazard[subset.Lambda0]
#'     hazard0.tindex <- Lambda0$hazard[subset.Lambda0]
#'     se.betaHazard.tindex <- matrix(NA, nrow = n.newdata, ncol = n.times)
#'   }
#'   if("cumHazard" %in% type || "survival" %in% type){
#'     se.cumHazard0.tindex <- Lambda0$se.cumHazard[subset.Lambda0]
#'     cumHazard0.tindex <- Lambda0$cumHazard[subset.Lambda0]
#'     se.betaCumHazard.tindex <- matrix(NA, nrow = n.newdata, ncol = n.times)
#'   }
#'   
#'   ##
#'   for(indexT in 1:n.times){
#'     
#'     if("hazard" %in% type){
#'       if(is.null(object$var)){
#'         se.betaHazard.tindex[,indexT] <- 0
#'       }else{
#'         Xbar_loop <- Lambda0$Xbar[subset.Lambda0[indexT],,drop = FALSE]
#'         hazard0_loop <- hazard0.tindex[indexT]
#'         
#'         hazardX_loop <- sweep(newdata*hazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(Xbar_loop))
#'         se.betaHazard.tindex[,indexT] <- rowSums(hazardX_loop %*% object$var * hazardX_loop)
#'       }
#'     }
#'     
#'     if("cumHazard" %in% type || "survival" %in% type){
#'       if(is.null(object$var)){
#'         se.betaCumHazard.tindex[,indexT] <- 0
#'       }else{
#'         XbarCumSum_loop <- Lambda0$XbarCumSum[subset.Lambda0[indexT],,drop = FALSE]
#'         cumHazard0_loop <- cumHazard0.tindex[indexT]
#'         cumHazardX_loop <- sweep(newdata*cumHazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(XbarCumSum_loop))
#'         se.betaCumHazard.tindex[,indexT] <- rowSums(cumHazardX_loop %*% object$var * cumHazardX_loop)  
#'       }
#'       
#'     }
#'     
#'   }
#'   
#'   ##
#'   out <- list()
#'   if("hazard" %in% type){
#'     tempo <- sqrt(sweep(se.betaHazard.tindex, MARGIN = 2, FUN = "+", STATS = se.hazard0.tindex))
#'     out$hazard.se <- sweep(tempo, MARGIN = 1, FUN = "*", STATS = eXb) 
#'   }
#'   if("cumHazard" %in% type){
#'     tempo <- sqrt(sweep(se.betaCumHazard.tindex, MARGIN = 2, FUN = "+", STATS = se.cumHazard0.tindex))
#'     out$cumHazard.se <- sweep(tempo, MARGIN = 1, FUN = "*", STATS = eXb) 
#'   }
#'   if("survival" %in% type){
#'     out$survival.se <- out$cumHazard.se*exp(-survival)
#'   }
#'   
#'   ## export
#'   return(out)
#' }
