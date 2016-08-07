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
#' @param keep.lastEventTime Logical. If \code{TRUE} add the time at which the last event occured for each strata to the output. 
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. Experimental !!
#' @details Not working with time varying predictor variables or
#'     delayed entry.
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
                       keep.strata = FALSE,
                       keep.times = FALSE,
                       keep.lastEventTime = FALSE,
                       se = FALSE){ 
  
  if(!is.null(newdata)){n.newdata <- NROW(newdata)}
  
  #### extract elements from objects ####
  xterms <- delete.response(object$terms)
  xvars <- attr(xterms,"term.labels")
  #### cph object
  if ("cph" %in% class(object)){
    nPatients <- sum(object$n)
    if(is.null(object$y)){
      stop("Argument \'y\' must be set to TRUE in cph \n")
    }
    strataspecials <- attr(xterms,"specials")$strat
    stratavars <- xvars[strataspecials]
    is.strata <- length(strataspecials)>0
    if(is.strata){ ## cph:strata for estimation of the baseline hazard
      if (length(xvars)>length(strataspecials)) 
        sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
      else 
        sterms <- xterms
      stratavars <- xvars[strataspecials]
      strataF <- object$Strata
    }else{
      strataF <- factor("1")
    }
    if(se){ ## cph:design matrix for standard error
      if(length(object$Design$mmcolnames)==0){
        #stop("Cannot compute standard errors when there is no confounder \n")
        modeldata <- matrix(0, nrow = nPatients, ncol = 1)
      }else{
        modeldata <- as.matrix(model.frame(object)[,object$Design$mmcolnames,drop=FALSE])
        modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$mean)
      }
    }else{
      modeldata <- matrix(0)
    }
    
  } else if ("coxph" %in% class(object)){ #### coxph object
    
    nPatients <- object$n
    strataspecials <- attr(xterms,"specials")$strata
    stratavars <- xvars[strataspecials]
    is.strata <- length(strataspecials)>0
    if(is.strata){ ## cph:strata for estimation of the baseline hazard
      if (length(xvars)>length(strataspecials)) 
        sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
      else 
        sterms <- xterms
      stratalevels <- object$xlevels[stratavars]
      strataF <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
    }else{
      strataF <- factor("1")
    }
    if(se){ ## cph:design matrix for standard error
      if(length(object$means)==0){
        stop("Cannot compute standard errors when there is no confounder \n")
        # modeldata <- matrix(0, nrow = nPatients, ncol = 1)
      }else{
        modeldata <- as.matrix(model.frame(object)[,names(object$means),drop = FALSE])
        modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$means)
      }
    }else{
      modeldata <- matrix(0)
    }
  } else {
    stop("Only implemented for \"coxph\" and \"cph\" objects \n")
  }
  ## checks
  if(object$method == "exact"){
    stop("Prediction with exact correction for ties is not implemented \n")
  }
  if(!missing(times) && any(is.na(times))){
    stop("NA values in argument \'times\' \n")
  }
  #   if(se && object$method != "efron"){
  #     stop("standard errors only implemented for object$method = efron \n")
  #   }
  if(!is.null(object$weights)){
    stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if(any(type %in% c("hazard","cumHazard","survival","eXb","newstrata") == FALSE)){
    stop("type can only be \"hazard\", \"cumHazard\" or/and \"survival\" \n") # eXb and newstrata are only for internal use 
  }
  
  #### baseline hazard ####
  levelsStrata <- levels(strataF)
  nStrata <- length(levelsStrata)
  ytimes <- object$y[,"time"]
  status <- object$y[,"status"]
  nVar <- ncol(modeldata)
  if(is.strata){ etimes.max <- tapply(ytimes, strataF, max) }else{ etimes.max <- max(ytimes) } # last event time
  
  Lambda0 <- baseHaz_cpp(alltimes = ytimes,
                         status = status,
                         eXb = if(centered == FALSE){exp(object$linear.predictors + sum(object$means*stats::coef(object)))}else{exp(object$linear.predictors)},
                         strata = as.numeric(strataF) - 1,
                         se = se,
                         data = modeldata,
                         nVar = nVar,
                         nPatients = nPatients,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = if(missing(times)){numeric(0)}else{sort(times)},
                         cause = 1,
                         Efron = (object$method == "efron"))
  
  if (is.strata == TRUE){ ## rename the strata value with the correct levels
    Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = levelsStrata)
  }
  
  #### linear predictor and strata for the new data ####
  if(!is.null(newdata)){ 
    if ("cph" %in% class(object)){
      if(length(xvars) > length(stratavars)){
        eXb <- exp(stats::predict(object, newdata, type = "lp"))
      }else{ 
        eXb <- rep(1, n.newdata) 
      }
      if(is.strata){
        tmp <- model.frame(sterms,newdata)
        colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
        tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
        newstrata <- apply(tmp,1,paste,collapse=".")
        newstrata <- factor(newstrata, levels = levelsStrata) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
      }
    } else if ("coxph" %in% class(object)){
      if(length(xvars) == length(stratavars)){ 
        eXb <- rep(1, n.newdata)
      } else if(is.strata){
        eXb <- exp(rowSums(stats::predict(object, newdata = newdata, type = "terms")))
      }else { 
        eXb <- exp(stats::predict(object, newdata, type = "lp"))
      }
      if(is.strata){
        newstrata <- prodlim::model.design(sterms,data=newdata,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
        newstrata <- factor(newstrata, levels = levelsStrata) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
      }
    }
    if(is.strata){
      allStrata <- unique(newstrata)
      if (any(allStrata %in% levelsStrata == FALSE)){
        stop("unknown strata: ",paste(unique(allStrata[allStrata %in% levelsStrata == FALSE]), collapse = " | "),"\n")
      }
    }
  }
  
  #### compute hazard and survival #### 
  if (is.null(newdata) || "eXb" %in% type || "newstrata" %in% type){  ## on the training dataset
    
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
    
    if ("eXb" %in% type && !is.null(newdata)){  # for internal use (called by predict.CauseSpecificCox), not documented
      Lambda0$eXb <- eXb
    }
    
    if ("newstrata" %in% type && !is.null(newdata)){ # for internal use (called by predict.CauseSpecificCox), not documented
      if (is.strata){Lambda0$newstrata = newstrata}else{Lambda0$newstrata = rep(1,n.newdata)}
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
    
    return(Lambda0)
    
  } else { ## on a new dataset
    out <- list()
    
    if(missing(times) || length(times)==0){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    n.times <- length(times)
    
    ## subject specific hazard
    if (is.strata==FALSE){
      if ("hazard" %in% type){out$hazard <- (eXb %o% Lambda0$hazard)}
      if ("cumHazard" %in% type || "survival" %in% type){
        cumHazard <- eXb %o% Lambda0$cumHazard
        if ("cumHazard" %in% type){out$cumHazard <- cumHazard}
        if ("survival" %in% type){out$survival <- exp(-cumHazard)}
      }
      
      if(se){ 
        out <- c(out,
                 seCox(object, newdata, times, type, Lambda0, survival = out$survival, eXb = eXb, stratavars = NULL)
        )
      }
      
    }else{ 
      
      ## initialization
      if ("hazard" %in% type){
        out$hazard <- matrix(0, nrow = n.newdata, ncol = n.times)
        if(se){out$hazard.se <- matrix(0, nrow = n.newdata, ncol = n.times)}
      }
      if ("cumHazard" %in% type){
        out$cumHazard <- matrix(NA, nrow = n.newdata, ncol = n.times)
        if(se){out$cumHazard.se <- matrix(0, nrow = n.newdata, ncol = n.times)}
      }
      if ("survival" %in% type){
        out$survival <- matrix(NA, nrow = n.newdata, ncol = n.times)
        if(se){out$survival.se <- matrix(0, nrow = n.newdata, ncol = n.times)}
      }
      
      ## loop across strata
      allStrata <- unique(newstrata)
      
      for(S in allStrata){
        id.S <- Lambda0$strata==S
        newid.S <- newstrata==S
        
        if ("hazard" %in% type){out$hazard[newid.S,] <- eXb[newid.S] %o% Lambda0$hazard[id.S]}
        if ("cumHazard" %in% type || "survival" %in% type){
          cumHazard.S <-  eXb[newid.S] %o% Lambda0$cumHazard[id.S]
          if ("cumHazard" %in% type){out$cumHazard[newid.S,] <- cumHazard.S}
          if ("survival" %in% type){out$survival[newid.S,] <- exp(-cumHazard.S)}
        }
        
        if(se){
          outSE <- seCox(object, newdata = newdata[newid.S,,drop = FALSE], 
                         times, type, 
                         Lambda0 = Lambda0, 
                         survival = out$survival[newid.S,], 
                         eXb = eXb[newid.S],
                         stratavars = stratavars, 
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
    if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=newstrata))
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
#' @param stratavars the name of the strata variables
#' @param subset.Lambda0 an index giving the subset of observations to be used in Lambda0
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
seCox <- function(object, newdata, times, type, Lambda0, eXb, survival, stratavars, subset.Lambda0 = 1:length(Lambda0$time)){
  
  ## prepare dataset
  if("cph" %in% class(object)){
    formulaObj <- object$sformula
    names.Xcenter <- object$Design$mmcolnames
    Xcenter <- object$mean
  }else if("coxph" %in% class(object)){
    formulaObj <- object$formula
    names.Xcenter <- names(object$means)
    Xcenter <- object$means
  }  
  
  terms.X <- delete.response(terms(formulaObj))
  if(!is.null(stratavars)){terms.X <- drop.terms(terms.X, which(attr(terms.X,"term.labels") %in% stratavars))}
  formulaX <- formula(terms.X)
  newdata.design <- stats::model.matrix(formulaX,newdata)[,names.Xcenter,drop=FALSE]
  newdata.design <- sweep(newdata.design, FUN = "-", MARGIN = 2, STATS = Xcenter)
  
  ##
  n.times <- length(times)
  n.newdata <- NROW(newdata)
  etimes <- Lambda0$time[subset.Lambda0]
  
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
      Xbar_loop <- Lambda0$Xbar[subset.Lambda0[indexT],,drop = FALSE]
      hazard0_loop <- hazard0.tindex[indexT]
      hazardX_loop <- sweep(newdata.design*hazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(Xbar_loop))
      se.betaHazard.tindex[,indexT] <- rowSums(hazardX_loop %*% object$var * hazardX_loop)
    }
    if("cumHazard" %in% type || "survival" %in% type){
      XbarCumSum_loop <- Lambda0$XbarCumSum[subset.Lambda0[indexT],,drop = FALSE]
      cumHazard0_loop <- cumHazard0.tindex[indexT]
      cumHazardX_loop <- sweep(newdata.design*cumHazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(XbarCumSum_loop))
      se.betaCumHazard.tindex[,indexT] <- rowSums(cumHazardX_loop %*% object$var * cumHazardX_loop)
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
