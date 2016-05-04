#' @title Predicting survival, hazard or cumulative hazard from Cox regression model 
#' 
#' @param object The fitted Cox regression model object either obtained with \code{coxph} (survival package) or \code{cph} (rms package).
#' @param newdata A data frame or table containing the values of the predictor variables 
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param type should the hazard or the cumulative hazard be returned
#' @param format should the result be returned as a matrix or a data.table
#' @param keep.strata If \code{TRUE} return an additional column containing the strata level. 
#' @param colnames Should the columns be named using the time to which they correspond. 
#' 
#' @details 
#' Note: for Cox regression models with time varying covariates this function may not work, and
#'       it does not make sense to ask for survival or cumulative hazard predictions.
#' 
#' @return A data table or a matrix containing the predictions for each subject (in rows)
#'         and each time (in columns) and the strata (if requested).
#' 
#' @examples 
#' 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d)
#' # table(duplicated(d$time))
#' ttt <- sample(x = unique(sort(d$time)), size = 10) 
#'
#' res1 <- predictSurv(fit, newdata = d, times = 10)
#' res2 <- predictSurv(fit, newdata = d, times = ttt)
#' res3 <- predictSurv(fit, newdata = d, times = 10, type = "cumHazard")
#' res4 <- predictSurv(fit, newdata = d, times = ttt, type = "cumHazard")
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- predictSurv(fitS, newdata = d, times = 10)
#' res2S <- predictSurv(fitS, newdata = d, times = ttt)
#' res3S <- predictSurv(fitS, newdata = d, times = 10, type = "cumHazard")
#' res4S <- predictSurv(fitS, newdata = d, times = ttt, type = "cumHazard")
#' @export
predictSurv <- function(object,
                        newdata,
                        times,
                        type = "Survival",
                        format = "data.table", 
                        keep.strata = FALSE,
                        colnames = FALSE) {
  
  match.arg(format, choices = c("matrix", "data.table"), several.ok = TRUE)
  nPatient <- NROW(newdata)  
  nTimes <- length(times)
  
  ## strata?
  if ("cph" %in% class(object))
    strataspecials <- attr(object$terms,"specials")$strat
  else
    strataspecials <- attr(object$terms,"specials")$strata 
  no.strata <- is.null(strataspecials)
  
  ## 
  BaseVar <- attr(object$terms, "term.labels")
  
  if(no.strata){ 
    
    stratalabels <- NULL
    strataVar <- NULL
    
  } else {
    
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    BaseVar <- BaseVar[ BaseVar %in% stratalabels == FALSE]
    
    if  ("cph" %in% class(object))
      strataVar <- names(prodlim:::parseSpecialNames(stratalabels, special = "strat"))
    else
      strataVar <- names(prodlim:::parseSpecialNames(stratalabels, special = "strata"))
  }
  
  ## linear predictor
  if  ("cph" %in% class(object)){
    if(length(BaseVar) > 0){
      Xb <- predict(object, newdata, type = "lp")}
    else{ 
      Xb <- rep(0, nrow(newdata)) 
    }
  }
  else{
    if(length(BaseVar) == 0){ 
      Xb <- rep(0, nrow(newdata))
    } else if(no.strata){ 
      Xb <- survival:::predict.coxph(object, newdata, type = "lp") 
    }else { 
      Xb <- rowSums(survival:::predict.coxph(object, newdata = newdata, type = "terms")) 
    }}
  
  ## baseline hazard
  Lambda0 = baseHazRR(object, centered = TRUE, maxtime =  max(times))

  ## preparation
  match.arg(type, choices = c("Survival", "cumHazard", "hazard"), several.ok = TRUE)
  n.type <- length(type)
  
  resPred <- list()
  if("Survival" %in% type){resPred$Survival <- matrix(NA, nrow = nPatient, ncol = nTimes)}else{resPred$Survival <- matrix(ncol=0,nrow=0)}
  if("cumHazard" %in% type){resPred$cumHazard <- matrix(NA, nrow = nPatient, ncol = nTimes)}else{resPred$cumHazard <- matrix(ncol=0,nrow=0)}
  if("hazard" %in% type){resPred$hazard <- matrix(0, nrow = nPatient, ncol = nTimes)}else{resPred$hazard <- matrix(ncol=0,nrow=0)}
  
  
  #### main
  if(no.strata){
    
      index.eventTimes <- intersect(which(times <= Lambda0[.N,time]),
                                    which(times >= 0))
      
      if(length(index.eventTimes)>0){
      
        if ("hazard" %in% type){
          resPred[["hazard"]][, which(times[index.eventTimes] %in% Lambda0$time)] <- Lambda0[which(Lambda0$time %in% times[index.eventTimes]), exp(Xb) %o% hazard]
        }
        if ("cumHazard" %in% type){
          jump.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times[index.eventTimes])
          resPred[["cumHazard"]][,index.eventTimes] <- exp(Xb) %o% c(0,Lambda0$cumHazard)[jump.index+1]
        }
        if ("Survival" %in% type){
          if("cumHazard" %in% type){
            resPred[["Survival"]] <- exp(- resPred[["cumHazard"]])
          }else{
            jump.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times[index.eventTimes])
            resPred[["Survival"]][,index.eventTimes] <- exp( - exp(Xb) %o% c(0,Lambda0$cumHazard)[jump.index+1] )
          }
        }
        
      }
     
       
#     resPred <- predictSurvStrata_cpp(cumHazard = Lambda0$cumHazard, 
#                                      hazard = Lambda0$hazard, 
#                                      eventtimes = Lambda0$time, 
#                                      times = times,
#                                      Xb = Xb,
#                                      originStrata = -1, 
#                                      newStrata = -1,
#                                      nStrata = 1,
#                                      f = prodlim::sindex,
#                                      returnSurvival = "Survival" %in% type, 
#                                      returnCumHazard = "cumHazard" %in% type, 
#                                      returnHazard = "hazard" %in% type)
  
    
  } else{ ## strata
     
    if(!is.null(attr(Xb,"strata"))){
      newStrata <- attr(Xb,"strata")
      newLevels <- levels(newStrata)
    } else {
      strataFormula <- as.formula(paste0("~ 0 + ",paste(strataVar,collapse = " + ")))
      newStrata <- interaction(model.frame(strataFormula,newdata), drop = TRUE, sep = ".")
      tempoLevels <- levels(newStrata)
      newLevels <- unlist(lapply(strsplit(tempoLevels, split = ".", fixed = TRUE), function(x){paste(paste(strataVar,x,sep = "="), collapse = ".")}))
    }
  
    originLevels <- levels(Lambda0$strata)
    
    if(any(newLevels %in% originLevels == FALSE)){
      stop("predictSurv_internal: unknown strata ",paste(newLevels[newLevels %in% originLevels == FALSE], collapse = " ")," \n",
           "existing strata: \"",paste(originLevels, collapse = "\" \""),"\"\n")
    }
    
    Lambda0$strata <- as.numeric(factor(Lambda0$strata, levels = c(newLevels, setdiff(originLevels,newLevels)))) - 1 
    newStrata <- as.numeric(newStrata) - 1
    
    ## compute hazard
      for(iterS in 0:(length(newLevels)-1)){
        Lambda0_Strata <- Lambda0[strata == iterS,]
        index.newPat <- which(newStrata == iterS)
        index.eventTimes <- intersect(which(times <= Lambda0_Strata[.N,time]),
                                      which(times >= 0))
        
        if(length(index.eventTimes) == 0){ next }
  
        if ("hazard" %in% type){
          resPred[["hazard"]][index.newPat, which(times[index.eventTimes] %in% Lambda0_Strata$time)] <- Lambda0_Strata[which(Lambda0_Strata$time %in% times[index.eventTimes]), exp(Xb[index.newPat]) %o% hazard]
        }
        if ("cumHazard" %in% type){
          jump.index <- prodlim::sindex(jump.times = Lambda0_Strata$time, eval.times = times[index.eventTimes])
          resPred[["cumHazard"]][index.newPat,index.eventTimes] <- exp(Xb[index.newPat]) %o% c(0,Lambda0_Strata$cumHazard)[jump.index+1]
        }
        if ("Survival" %in% type){
          if("cumHazard" %in% type){
            resPred[["Survival"]][index.newPat,] <- exp(- resPred[["cumHazard"]][index.newPat,])
          }else{
            jump.index <- prodlim::sindex(jump.times = Lambda0_Strata$time, eval.times = times[index.eventTimes])
            resPred[["Survival"]][index.newPat,index.eventTimes] <- exp(- exp(Xb[index.newPat]) %o% c(0,Lambda0_Strata$cumHazard)[jump.index+1] )
          }
        }
        
      }
      
#      resPred <- predictSurvStrata_cpp(cumHazard = Lambda0$cumHazard, 
#                                        hazard = Lambda0$hazard, 
#                                        eventtimes = Lambda0$time, 
#                                        times = times,
#                                        Xb = Xb,
#                                        originStrata = Lambda0$strata, 
#                                        newStrata = newStrata,
#                                        nStrata = length(newLevels),
#                                        f = prodlim::sindex,
#                                        returnSurvival = "Survival" %in% type, 
#                                        returnCumHazard = "cumHazard" %in% type, 
#                                        returnHazard = "hazard" %in% type)
  
  }
  
  ## export
  if("Survival" %in% type){
    if(colnames){colnames(resPred$Survival) <- as.character(times)}
    if(format == "data.table"){resPred$Survival <- as.data.table(resPred$Survival)}
    
    if(length(type)==1){
      if(format == "data.table" && keep.strata == TRUE && length(strataVar) > 0){
        return(data.table(strata = strata, resPred$Survival))
      }else{
        return(resPred$Survival)
      }
    }
  }
  if("cumHazard" %in% type){
    if(colnames){colnames(resPred$cumHazard) <- as.character(times)}
    if(format == "data.table"){resPred$cumHazard <- as.data.table(resPred$cumHazard)}
    
    if(length(type)==1){
      if(format == "data.table" && keep.strata == TRUE && length(strataVar) > 0){
        return(data.table(strata = strata, resPred$cumHazard))
      }else{
        return(resPred$cumHazard)
      }
    }
  }
  if("hazard" %in% type){
    if(colnames){colnames(resPred$hazard) <- as.character(times)}
    if(format == "data.table"){resPred$hazard <- as.data.table(resPred$hazard)}
    
    if(length(type)==1){
      if(format == "data.table" && keep.strata == TRUE && length(strataVar) > 0){
        return(data.table(strata = strata, resPred$hazard))
      }else{
        return(resPred$hazard)
      }
    }
  }
  
  return(resPred)
  
}




