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
  
  if(!is.null(newdata)){n.newdata <- NROW(newdata)}
  
  #### extract elements from objects ####
  resInfo <- getCoxInfo(object, design = se)
  is.strata <- resInfo$is.strata
 
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
  if(any(type %in% c("hazard","cumHazard","survival","eXb","newstrata") == FALSE)){
    stop("type can only be \"hazard\", \"cumHazard\" or/and \"survival\" \n") # eXb and newstrata are only for internal use 
  }
  # if(se == TRUE && ncol(resInfo$modeldata) == 0){
  #   stop("cannot compute the standard error when there are not covariates \n")
  # }
  
  #### baseline hazard ####
  levelsStrata <- levels(resInfo$strataF)
  nStrata <- length(levelsStrata)
  ytimes <- object$y[,"time"]
  status <- object$y[,"status"]
  nVar <- ncol(resInfo$modeldata)
  if(is.strata){ etimes.max <- tapply(ytimes, resInfo$strataF, max) }else{ etimes.max <- max(ytimes) } # last event time
  
  Lambda0 <- baseHaz_cpp(alltimes = ytimes,
                         status = status,
                         eXb = if(centered == FALSE){exp(object$linear.predictors + sum(object$means*stats::coef(object)))}else{exp(object$linear.predictors)},
                         strata = as.numeric(resInfo$strataF) - 1,
                         se = se,
                         data = resInfo$modeldata,
                         nVar = nVar,
                         nPatients = resInfo$nPatients,
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

      eXb <- exp(lpCox(object, data = newdata, 
                   xvars = resInfo$xvars, stratavars = resInfo$stratavars))
              
      if(is.strata){
        newstrata <- defineStrata(object, data = newdata, 
                                  sterms = resInfo$sterms, 
                                  stratavars = resInfo$stratavars, 
                                  levelsStrata = levelsStrata, 
                                  stratalevels = resInfo$stratalevels)
        
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
    
    return(do.call(paste0("as.",format), list(Lambda0)))
    
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
                         stratavars = resInfo$stratavars, 
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
  newdata.design <- centerData(object, data = newdata, stratavars = stratavars)
 
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
      if(is.null(object$var)){
        se.betaHazard.tindex[,indexT] <- 0
      }else{
        Xbar_loop <- Lambda0$Xbar[subset.Lambda0[indexT],,drop = FALSE]
        hazard0_loop <- hazard0.tindex[indexT]
        hazardX_loop <- sweep(newdata.design*hazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(Xbar_loop))
        se.betaHazard.tindex[,indexT] <- rowSums(hazardX_loop %*% object$var * hazardX_loop)
      }
    }
    
    if("cumHazard" %in% type || "survival" %in% type){
      if(is.null(object$var)){
        se.betaCumHazard.tindex[,indexT] <- 0
      }else{
        XbarCumSum_loop <- Lambda0$XbarCumSum[subset.Lambda0[indexT],,drop = FALSE]
        cumHazard0_loop <- cumHazard0.tindex[indexT]
        cumHazardX_loop <- sweep(newdata.design*cumHazard0_loop, MARGIN = 2, FUN = "-",  STATS = as.double(XbarCumSum_loop))
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


#' @title Extract information from a Cox model
#' @description Extract information from a Cox model that are necessary for computing the baseline hazard, performing predictions or estimating the influence function associated to the Cox model
#' @name getCoxInfo
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param design should the design matrix be returned? 
#'
#' @return A named list containing the following elements:
#' \itemize{
#'  \item{"nPatients"}{the number of observations}
#'  \item{"xvars"}{the name of all the regressors (including those used to form the strata)}
#'  \item{"modeldata"}{the design matrix for the regressors}
#'  \item{"stratavars"}{the name of the variables used to define the strata}
#'  \item{"is.strata"}{is there any strata?}
#'  \item{"sterms"}{terms corresponding to the strata variables}
#'  \item{"strataF"}{a vector contain the strata factor for each observation}
#'  \item{"stratalevels"}{a named list containing for each variable used to form the strata all its possible levels}
#' }
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' getCoxInfo(mCox, design = FALSE) 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' getCoxInfo(mCoxS, design = FALSE) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' getCoxInfo(mCox, design = FALSE) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' getCoxInfo(mCoxS, design = FALSE) 
#' }

#' @rdname getCoxInfo
getCoxInfo <- function(object, design) UseMethod("getCoxInfo")

#' @rdname getCoxInfo
getCoxInfo.cph <- function(object, design){
  
  xterms <- delete.response(object$terms)
  xvars <- attr(xterms,"term.labels")
  
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
    strataF <- object$strata
  }else{
    sterms <- NULL
    stratalevels <- NULL
    strataF <- factor("1")
  }
  if(design){ ## cph:design matrix for standard error
    modeldata <- centerData(object, data = model.frame(object), stratavars = stratavars)
  }else{
    modeldata <- matrix(nrow = 0, ncol = 0)
  }
  
  return(list(nPatients = nPatients,
              xvars = xvars,
              modeldata = modeldata,
              stratavars = stratavars, 
              is.strata = is.strata,
              strataF = strataF,
              stratalevels = NULL,
              sterms = sterms))
}

#' @rdname getCoxInfo
getCoxInfo.coxph <- function(object, design){
  
  xterms <- delete.response(object$terms)
  xvars <- attr(xterms,"term.labels")
  
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
    sterms <- NULL
    stratalevels <- NULL
  }
  
  
  if(design){ ## coxph:design matrix for standard error
    modeldata <- centerData(object, data = model.frame(object), stratavars = stratavars)
  }else{
    modeldata <- matrix(nrow = 0, ncol = 0)
  }
  
  return(list(nPatients = nPatients,
              xvars = xvars,
              modeldata = modeldata,
              stratavars = stratavars, 
              is.strata = is.strata,
              strataF = strataF,
              stratalevels = stratalevels,
              sterms = sterms))
}


#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname lpCox 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param xvars the name of all the regressors (including those used to form the strata)
#' @param stratavars the variables used to define the strata
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL) 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' lpCox(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "X1") 
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' lpCox(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2")) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' lpCox(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "X1") 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' lpCox(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2")) 
#' }

#' @rdname lpCox
lpCox <- function(object, data, xvars, stratavars) UseMethod("lpCox")

#' @rdname lpCox
lpCox.cph <- function(object, data, xvars, stratavars){
  
  if(length(xvars) > length(stratavars)){
    Xb <- stats::predict(object, data, type = "lp")
  }else{ 
    Xb <- rep(0, NROW(data)) 
  }
  
  return(Xb)
}
  
#' @rdname lpCox
lpCox.coxph <- function(object, data, xvars, stratavars){

  if(length(xvars) == length(stratavars)){ 
    Xb <- rep(0, NROW(data))
  } else if(length(stratavars)>0){
    Xb <- rowSums(stats::predict(object, newdata = data, type = "terms"))
  }else { 
    Xb <- stats::predict(object, data, type = "lp")
  }
  
  return(Xb)
}

#' @title Define the strata for a new dataset
#' @description Define the strata in a dataset to match those of a stratified Cox model
#' @name defineStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param sterms terms in the formula corresponding to the strata variables
#' @param stratavars the name of the variables used to define the strata
#' @param levelsStrata the strata levels that have been used to fit the Cox model
#' @param stratalevels a named list containing for each variable used to form the strata all its possible levels
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' info <- getCoxInfo(mCoxS, design = FALSE)
#' defineStrata(mCoxS, data = d, sterms = info$sterms, stratavars = info$stratavars, 
#'              levelsStrata = levels(info$strataF), stratalevels = info$stratalevels) 
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d, y = TRUE)
#' info <- getCoxInfo(mCoxS, design = FALSE)
#' defineStrata(mCoxS, data = d, sterms = info$sterms, stratavars = info$stratavars, 
#'              levelsStrata = levels(info$strataF), stratalevels = info$stratalevels) 
#' }
defineStrata <- function(object, data, sterms, stratavars, levelsStrata, ...) UseMethod("defineStrata")


#' @rdname defineStrata
defineStrata.cph <- function(object, data, sterms, stratavars, levelsStrata, ...){

    if(length(stratavars)>0){
    tmp <- model.frame(sterms,data)
    colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
    tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
    newstrata <- apply(tmp,1,paste,collapse=".")
    newstrata <- factor(newstrata, levels = levelsStrata) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
  }else{
    newstrata <- NULL
  } 
  
  return(newstrata)
}
  
#' @rdname defineStrata
defineStrata.coxph <- function(object, data, sterms, stratavars, levelsStrata, stratalevels){
  
  if(length(stratavars)>0){
    newstrata <- prodlim::model.design(sterms,data=data,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
    newstrata <- factor(newstrata, levels = levelsStrata) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
  }else{
    newstrata <- NULL
  } 
  
  return(newstrata)
}
  
   

#' @title Center a dataset
#' @description Center a dataset around the centering values
#' @name centerData
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ X1, data = d) 
#' centerData(mCoxS, stratavars = NULL)
#' mCoxS <- coxph(Surv(time, event) ~ X1+X2, data = d) 
#' centerData(mCoxS, stratavars = NULL)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d) 
#' centerData(mCoxS, stratavars = "strata(X1)")
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d) 
#' centerData(mCoxS, stratavars = c("strata(X1)","strata(X2)"))
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ X1, data = d, y = TRUE) 
#' centerData(mCoxS, stratavars = NULL)
#' mCoxS <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE) 
#' centerData(mCoxS, stratavars = NULL)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d, y = TRUE) 
#' centerData(mCoxS, stratavars = "strat(X1)")
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE) 
#' centerData(mCoxS, stratavars = c("strat(X1)","strat(X2)"))
#' }
centerData <- function(object, data, stratavars, center, ...) UseMethod("centerData")

#' @rdname centerData
centerData.cph <- function(object, data = model.frame(object), stratavars, center = TRUE, ...){

  if(length(object$Design$mmcolnames)>0){
    
    terms.X <- delete.response(terms(object$sformula))
    if(length(stratavars)>0){terms.X <- drop.terms(terms.X, which(attr(terms.X,"term.labels") %in% stratavars))}
    formulaX <- formula(terms.X)
    modeldata <- stats::model.matrix(formulaX,data)[,object$Design$mmcolnames,drop=FALSE]
    
    if(center){
      modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$means)
    }
    
  }else{
    
    modeldata <- matrix(0, ncol = 1, nrow = NROW(data))
    
  }
  
  return(modeldata)
}

#' @rdname centerData
centerData.coxph <- function(object, data = model.frame(object), stratavars, center = TRUE, ....){
  
  if(length(object$means)>0){
    
    terms.X <- delete.response(terms(object$formula))
    if(length(stratavars)>0){terms.X <- drop.terms(terms.X, which(attr(terms.X,"term.labels") %in% stratavars))}
    formulaX <- formula(terms.X)
    modeldata <- stats::model.matrix(formulaX,data)[,names(object$means),drop=FALSE]
   
    if(center){
      modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$means)
    }
    
  }else{
    
    modeldata <- matrix(0, ncol = 1, nrow = NROW(data))
    
  }
  
  return(modeldata)
}
