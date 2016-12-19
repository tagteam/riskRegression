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
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. Experimental !!
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
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#' 
#' predictCox(fit)
#' predictCox(fit, newdata=nd, times = 5)
#' cbind(survival::basehaz(fit),predictCox(fit,type="cumHazard"))
#' 
#' # one strata variable
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
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
#' fit2S <- coxph(Surv(time,status)~X1+strata(U)+strata(V)+X2,data=d, ties="breslow")
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
                       format = "data.frame"){ 
  
  #### extract elements from object ####
  # we need:          - the total number of observations, the status and eventtime for each observation
  #                   - the strata corresponding to each observation in the training set
  #                   - the name of each strata
  #                   - the value of the linear predictor for each observation in the training set
  # (optional for se) - the centered dataset for the variables involved the linear predictor
  infoVar <- CoxStrataVar(object)
  is.strata <- infoVar$is.strata
  
  object.n <- CoxN(object)
  object.status <- CoxStatus(object)
  object.time <- CoxEventtime(object)
  object.strata <- CoxStrata(object, stratavars = infoVar$stratavars)
  object.levelStrata <- levels(object.strata)
  object.eXb <- exp(CoxLP(object, data = NULL, center = centered))
  nVar <- length(infoVar$lpvars)
  
  if(se){
    object.LPdata <- as.matrix(CoxDesign(object, data = CoxData(object), 
                                         lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                         rm.intercept = TRUE, center = TRUE))
  }else{
    object.LPdata <- matrix(ncol = 0, nrow = 0)
  }
  
  #### Do we want to make prediction for a new dataset ? ####
  # if yes we need to define: - the linear predictor for the new dataset
  #                           - the strata corresponding to each observation in the new dataset
  # (optional for se)         - subset the dataset to the variable involved in the linear predictor and center these variables
  if(!is.null(newdata)){
    new.n <- NROW(newdata)
    newdata <- as.data.table(newdata)
     new.eXb <- exp(CoxLP(object, data = newdata, center = FALSE)) # must be the original name of the strata variables
    
    new.strata <- CoxStrata(object, data = newdata, 
                            sterms = infoVar$sterms, 
                            stratavars = infoVar$stratavars, 
                            levels = object.levelStrata, 
                            stratalevels = infoVar$stratalevels)
    
    new.levelStrata <- levels(new.strata)
    
    if(se){
      new.LPdata <- as.matrix(CoxDesign(object, data = newdata, 
                                        lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                        rm.intercept = TRUE, center = TRUE))
    }
    
  }
  
  ## checks
  if(object$method == "exact"){
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
  # if(se == TRUE && ncol(resInfo$modeldata) == 0){
  #   stop("cannot compute the standard error when there are not covariates \n")
  # }
  
  #### baseline hazard ####
  nStrata <- length(object.levelStrata)
  if(is.strata){etimes.max <- tapply(object.time, object.strata, max) }else{ etimes.max <- max(object.time) } # last event time
  
   Lambda0 <- baseHaz_cpp(alltimes = object.time,
                         status = object.status,
                         eXb = object.eXb,
                         strata = as.numeric(object.strata) - 1,
                         se = se,
                         data = object.LPdata,
                         nVar = max(1,nVar), # if no LP we still use one constant variable with value = 0
                         nPatients = object.n,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = if(missing(times)){numeric(0)}else{sort(times)},
                         cause = 1,
                         Efron = (object$method == "efron"))
   
  if (is.strata == TRUE){ ## rename the strata value with the correct levels
    Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
  }
  
  #### compute hazard and survival #### 
  if (is.null(newdata)){  ## on the training dataset
    
    Lambda0$Xbar <- NULL
    Lambda0$XbarCumSumRes <- NULL
    
    if ("hazard" %in% type == FALSE){ Lambda0$hazard <- NULL } 
    if ("hazard" %in% type == FALSE || se == FALSE){  Lambda0$se.hazard <- NULL }
    
    if ("cumHazard" %in% type == FALSE){ Lambda0$cumHazard <- NULL } 
    if ("cumHazard" %in% type == FALSE || se == FALSE){  Lambda0$se.cumHazard <- NULL }
    
    if ("survival" %in% type){ 
      Lambda0$survival = exp(-Lambda0$cumHazard)
      if(se){Lambda0$se.survival = sqrt( Lambda0$se.cumHazard*exp(-2*Lambda0$cumHazard) )} # delta method: Var(exp(b)) = exp(2b)var(b)
    } 
    
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
      
      if(se){ 
        out <- c(out,
                 seCox(object, newdata = new.LPdata, times = times, type = type, 
                       Lambda0 = Lambda0, eXb = new.eXb, survival = out$survival, 
                       subset.Lambda0 = 1:length(Lambda0$time))
        )
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
        
        if(se){
          outSE <- seCox(object, newdata = new.LPdata[newid.S,,drop = FALSE], 
                         times, type, 
                         Lambda0 = Lambda0, 
                         survival = out$survival[newid.S,], 
                         eXb = new.eXb[newid.S],
                         subset.Lambda0 = which(id.S))
          
          if ("hazard" %in% type){out$hazard.se[newid.S,] <- outSE$hazard.se}
          if ("cumHazard" %in% type){out$cumHazard.se[newid.S,] <- outSE$cumHazard.se}
          if ("survival" %in% type){out$survival.se[newid.S,] <- outSE$survival.se}
        }
      }
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
#' @param newdata  A data frame or table containing the values of the
#'     predictor variables defining subject specific predictions. Should have
#'     the same structure as the data set used to fit the \code{object}.
#' @param times Time points at which to evaluate the predictions. 
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param Lambda0 the baseline hazard estimate returned by BaseHazStrata_cpp.
#' @param eXb exponential (covariates * coefficients)
#' @param survival the predicted survival (only necessary when type contains survival)
#' @param subset.Lambda0 an index giving the subset of observations to be used in Lambda0
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
seCox <- function(object, newdata, times, type, Lambda0, eXb, survival, subset.Lambda0 = 1:length(Lambda0$time)){
  
  ## WARNING newdata must be centered
  n.times <- length(times)
  n.newdata <- NROW(newdata)
  etimes <- Lambda0$time[subset.Lambda0]
  newdata <- as.matrix(newdata)
  
  ##
  if("hazard" %in% type){
    se.hazard0.tindex <- Lambda0$se.hazard[subset.Lambda0]
    hazard0.tindex <- Lambda0$hazard[subset.Lambda0]
    se.betaHazard.tindex <- matrix(NA, nrow = n.newdata, ncol = n.times)
  }
  if("cumHazard" %in% type || "survival" %in% type){
    se.cumHazard0.tindex <- Lambda0$se.cumHazard[subset.Lambda0]
    cumHazard0.tindex <- Lambda0$cumHazard[subset.Lambda0]
    se.betaCumHazard.tindex <- matrix(NA, nrow = n.newdata, ncol = n.times)
  }
  
  ##
  for(indexT in 1:n.times){
    
    if("hazard" %in% type){
      if(is.null(object$var)){
        se.betaHazard.tindex[,indexT] <- 0
      }else{
        Xbar_loop <- Lambda0$Xbar[subset.Lambda0[indexT],,drop = FALSE]
        hazard0_loop <- hazard0.tindex[indexT]
        
        hazardX_loop <- sweep(newdata*hazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(Xbar_loop))
        se.betaHazard.tindex[,indexT] <- rowSums(hazardX_loop %*% object$var * hazardX_loop)
      }
    }
    
    if("cumHazard" %in% type || "survival" %in% type){
      if(is.null(object$var)){
        se.betaCumHazard.tindex[,indexT] <- 0
      }else{
        XbarCumSum_loop <- Lambda0$XbarCumSum[subset.Lambda0[indexT],,drop = FALSE]
        cumHazard0_loop <- cumHazard0.tindex[indexT]
        cumHazardX_loop <- sweep(newdata*cumHazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(XbarCumSum_loop))
        se.betaCumHazard.tindex[,indexT] <- rowSums(cumHazardX_loop %*% object$var * cumHazardX_loop)  
      }
      
    }
    
  }
  
  ##
  out <- list()
  if("hazard" %in% type){
    tempo <- sqrt(sweep(se.betaHazard.tindex, MARGIN = 2, FUN = "+", STATS = se.hazard0.tindex))
    out$hazard.se <- sweep(tempo, MARGIN = 1, FUN = "*", STATS = eXb) 
  }
  if("cumHazard" %in% type){
    tempo <- sqrt(sweep(se.betaCumHazard.tindex, MARGIN = 2, FUN = "+", STATS = se.cumHazard0.tindex))
    out$cumHazard.se <- sweep(tempo, MARGIN = 1, FUN = "*", STATS = eXb) 
  }
  if("survival" %in% type){
    out$survival.se <- out$cumHazard.se*exp(-survival)
  }
  
  ## export
  return(out)
}





