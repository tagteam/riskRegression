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
#' @param maxtime Baseline hazard will be computed for each event before maxtime
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}
#' @param keep.strata Logical. If \code{TRUE} add the (newdata) strata to the output. Only if there any. 
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output. 
#' @param keep.lastEventTime Logical. If \code{TRUE} add the time at which the last event occured for each strata to the output. 
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. 
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
#' predictCox(fit, maxtime = 5,keep.times=TRUE)
#' predictCox(fit, maxtime = 1,keep.times=TRUE)
#' 
#' # one strata variable
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' predictCox(fitS)
#' predictCox(fitS, maxtime = 5)
#' predictCox(fitS, maxtime = 5, newdata=nd, times = 1)
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
#' predictCox(fitS, maxtime = 5)
#' predictCox(fitS, maxtime = 5,newdata=nd, times = 3)
#' 
#' 
#' @export
predictCox <- function(object,
                       newdata=NULL,
                       times,
                       centered = TRUE,
                       maxtime = NULL,
                       type=c("hazard","cumHazard","survival"),
                       keep.strata = FALSE,
                       keep.times = FALSE,
                       keep.lastEventTime = FALSE,
                       se = FALSE){ 
  ## extract elements from objects
  xterms <- delete.response(object$terms)
  xvars <- attr(xterms,"term.labels")
  if ("cph" %in% class(object)){
    nPatients <- sum(object$n)
    if(is.null(object$y)){
      stop("Argument \'y\' must be set to TRUE in cph \n")
    }
    strataspecials <- attr(xterms,"specials")$strat
    stratavars <- xvars[strataspecials]
    is.strata <- length(strataspecials)>0
    if(is.strata){
      if (length(xvars)>length(strataspecials)) 
        sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
      else 
        sterms <- xterms
      stratavars <- xvars[strataspecials]
      strataF <- object$Strata
      ## stratalevels <- levels(factor(object$strata))
      ## names(stratalevels) <- stratavars
    }else{
      strataF <- factor("1")
    }
  } else if ("coxph" %in% class(object)){
    nPatients <- object$n
    strataspecials <- attr(xterms,"specials")$strata
    stratavars <- xvars[strataspecials]
    is.strata <- length(strataspecials)>0
    if(is.strata){
      if (length(xvars)>length(strataspecials)) 
        sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
      else 
        sterms <- xterms
      stratalevels <- object$xlevels[stratavars]
      strataF <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
    }else{
      strataF <- factor("1")
    }
  } else {
    stop("Only implemented for \"coxph\" and \"cph\" objects \n")
  }
  ## method
  if(object$method == "exact"){
    stop("Prediction with exact correction for ties is not implemented \n")
  }
  if(!missing(times) && any(is.na(times))){
    stop("NA values in argument \'times\' \n")
  }
  if(se && object$method != "efron"){
    stop("standard errors only implemented for object$method = efron \n")
  }
  ## main
  levelsStrata <- levels(strataF)
  nStrata <- length(levelsStrata)
  ytimes <- object$y[,"time"]
  status <- object$y[,"status"]
  
  if (is.null(maxtime)){ # avoid useless computation of the baseline hazard for times after the prediction time
    if(!is.null(newdata) && !missing(times)){
      etimes.min <- if(is.strata){tapply(ytimes, strataF, min)}else{min(ytimes)} # first event in each strata
      maxtime <- max(c(times,etimes.min)) + 1e-5
    }else{
      maxtime <- Inf
    }
  } 
  
  Lambda0 <- BaseHazStrata_cpp(alltimes = ytimes,
                               status = status,
                               Xb = if(centered == FALSE){object$linear.predictors + sum(object$means*stats::coef(object))}else{object$linear.predictors},
                               strata = as.numeric(strataF) - 1,
                               nPatients = nPatients,
                               nStrata = nStrata,
                               maxtime = maxtime,
                               cause = 1,
                               Efron = (object$method == "efron"),
                               se = se)
  
  if (is.strata == TRUE){ ## rename the strata value with the correct levels
    Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = levelsStrata)
  }
  
  out <- list()
  if (is.null(newdata)){ 
    
    if ("hazard" %in% type){ 
      out <- c(out, list(hazard = Lambda0$hazard))
      if(se){out <- c(out, list(se.hazard = Lambda0$se.hazard))}
    } 
    if ("cumHazard" %in% type){ 
      out <- c(out, list(cumHazard = Lambda0$cumHazard))
      if(se){out <- c(out, list(se.cumHazard = Lambda0$se.cumHazard))}
    } 
    if ("survival" %in% type){ 
      out <- c(out, list(survival =exp(-Lambda0$cumHazard)))
      if(se){out <- c(out,list(se.survival = NA))}
    } 
    if (keep.times==TRUE){
      out <- c(out,list(times=Lambda0$time))
    } 
    if(is.strata && keep.strata==TRUE){
      newstrata <- Lambda0$strata
    }
    
  } else {
    
    if(missing(times)){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    n.times <- length(times)
    
    ## linear predictor
    if  ("cph" %in% class(object)){
      if(length(xvars) > length(stratavars)){
        Xb <- stats::predict(object, newdata, type = "lp")}
      else{ 
        Xb <- rep(0, NROW(newdata)) 
      }
    } else{ ## coxph
      if(length(xvars) == length(stratavars)){ 
        Xb <- rep(0, NROW(newdata))
      } else if(is.strata){
        Xb <- rowSums(stats::predict(object, newdata = newdata, type = "terms")) 
      }else { 
        Xb <- stats::predict(object, newdata, type = "lp") 
      }
    }
    
    findInterval2 <- function(...){ # if all the prediction times are before the first event then return the index of the first event
      res <- setdiff(findInterval(...),0)
      if(length(res)==0){return(1)}else{res}
    }
    
    ## subject specific hazard
    if (is.strata==FALSE){
      
      ## remove useless event times for prediction, i.e. keep only event times that are the closest non-superior time for a prediction time
      keep.eventtime <- unique(sort(findInterval2(times, vec = Lambda0$time)))  # keep ascending time order
      Lambda0 <- lapply(Lambda0, function(x){x[keep.eventtime]})
      
      etimes <- Lambda0$time
      etimes.max <- max(ytimes) # last event time (cannot be max(Lambda0$time) because the censored observations have been removed from Lambda0)
      test.timeNA <- times>etimes.max
      
      if ("hazard" %in% type){
        hits <- times%in%etimes
        if (sum(hits)==0) { ## no time corresponds to an event
          out$hazard <- matrix(0,ncol=n.times,nrow=NROW(newdata))
        }else{  ## some times correspond to an event
          out$hazard <- matrix(0, nrow = NROW(newdata), ncol = n.times)
          out$hazard[,hits] <- exp(Xb) %o% Lambda0$hazard[match(times[hits], etimes)] #  match is needed here instead of %in% to handle non-increasing times.
        }
        if(any(test.timeNA)){out$hazard[,times>etimes.max] <- NA}
      }
      
      if ("cumHazard" %in% type || "survival" %in% type){
        cumHazard <- exp(Xb) %o% Lambda0$cumHazard
        
        tindex <- prodlim::sindex(jump.times=etimes,eval.times=times)
        if(any(test.timeNA)){tindex[times>etimes.max] <- NA}
      }
      
      if ("cumHazard" %in% type){
        out$cumHazard <- cbind(0,cumHazard)[,tindex+1, drop = FALSE]
      }
      
      if ("survival" %in% type){
        out$survival <- cbind(1,exp(-cumHazard))[,tindex+1, drop = FALSE]
      }
      
    }else{ 
      ## remove useless event times for prediction, i.e. keep only event times that are, within each strata, the closest non-superior time for a prediction time
      ## index[1] - 1  is 0 for the first strata then the index of the observation just before the second strata and so on
     keep.eventtime <- tapply(1:length(Lambda0$strata), Lambda0$strata, 
                               function(index){index[1] - 1 + findInterval2(times, vec = Lambda0$time[index])})
       keep.eventtime <- unique(sort(unlist(keep.eventtime))) # keep ascending time/strata order
      Lambda0 <- lapply(Lambda0, function(x){x[keep.eventtime]})
      
      etimes.max <- tapply(ytimes, strataF, max) # last event time (cannot be tapply(Lambda0$times, Lambda0$strata, tail, n = 1) because the censored observations have been removed from Lambda0)
     
      ## find strata for the new observations
      if ("cph" %in% class(object)){
        tmp <- model.frame(sterms,newdata)
        colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
        tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
        newstrata <- apply(tmp,1,paste,collapse=".")
        ## newstrata <- prodlim::model.design(sterms,data=newdata,specialsFactor=TRUE)$strat[[1]]                
      }else{
        newstrata <- prodlim::model.design(sterms,data=newdata,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
      }
      ## initialization
      if ("hazard" %in% type){out$hazard <- matrix(0, nrow = NROW(newdata), ncol = n.times)}
      if ("cumHazard" %in% type){out$cumHazard <- matrix(NA, nrow = NROW(newdata), ncol = n.times)}
      if ("survival" %in% type){out$survival <- matrix(NA, nrow = NROW(newdata), ncol = n.times)}
      
      ## loop across strata
      allStrata <- unique(newstrata)
      if (any(allStrata %in% levelsStrata == FALSE)){
        stop("unknown strata: ",paste(unique(allStrata[allStrata %in% levelsStrata == FALSE]), collapse = " | "),"\n")
      }
      for(S in allStrata){
        id.S <- Lambda0$strata==S
        newid.S <- newstrata==S
        etimes.S <- Lambda0$time[id.S]
        test.timeNA <- times>etimes.max[S]
        ## 
        if ("hazard" %in% type){
          hits <- times%in%etimes.S
          if (sum(hits)==0) {
            out$hazard[newid.S,] <- matrix(0,nrow=sum(newid.S),ncol=n.times)
          }else{
            out$hazard[newid.S,hits] <- exp(Xb[newid.S]) %o% Lambda0$hazard[id.S][match(times[hits], etimes.S)] #  match is needed here instead of %in% to handle non-increasing times.
          }
          if(any(test.timeNA)){out$hazard[newid.S,test.timeNA] <- NA}
        }
        
        if ("cumHazard" %in% type || "survival" %in% type){
          tindex.S <- prodlim::sindex(jump.times=etimes.S,eval.times=times)
          if(any(test.timeNA)){tindex.S[test.timeNA] <- NA}
          cumHazard.S <-  exp(Xb[newid.S]) %o% Lambda0$cumHazard[id.S]
        }
        if ("cumHazard" %in% type){
          out$cumHazard[newid.S,] <- cbind(0,cumHazard.S)[,tindex.S+1, drop = FALSE]
        }
        if ("survival" %in% type){
          out$survival[newid.S,] <- cbind(1,exp(-cumHazard.S))[,tindex.S+1, drop = FALSE]
        } 
      }
    }
    
    if (keep.times==TRUE) out <- c(out,list(times=times))
  }
  
  if( keep.lastEventTime==TRUE) out <- c(out,list(lastEventTime=etimes.max))
  if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=newstrata))
  return(out)
}

