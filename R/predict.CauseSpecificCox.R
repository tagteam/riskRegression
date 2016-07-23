#' Predicting hazard or cumulative hazard
#'
#' Apply formula to combine two or more Cox models into absolute risk
#' (cumulative incidence function)
#' @param object The fitted cause specific Cox model
#' @param newdata A data frame containing the values of the variables
#'     in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated
#'     hazards/survival
#' @param cause should the hazard or the cumulative hazard be returned
#' @param keep.strata If \code{TRUE} return an additional column
#'     containing the strata level.
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. Experimental !!
#' @param ... not used
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds
#'     tag@@biostat.ku.dk
#' @details Note: for Cox regression models with time varying
#'     covariates it does not make sense to use this function, because
#'     the predicted risk has to be a measurable function of the data
#'     available at the time origin.
#' 
#' @return A data table containing the predictions for each patient (in rows) and each time (in columns). 
#' 
#' @examples 
#' 
#' d <- SimCompRisk(1e2)
#' d$time <- round(d$time,1)
#' ttt <- sample(x = unique(sort(d$time)), size = 10) 
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow" )
#' 
#' predCSC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
#' 
#' #### cph function
#' CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", fitter = "cph")
#' 
#' predcph <- predict(CSC.cph, newdata = d, cause = 2, times = ttt)
#'
#' 
#' 
#' @export
predict.CauseSpecificCox <- function(object,newdata,times,cause,keep.strata = FALSE, se = FALSE, ...){
    survtype <- object$survtype
    if (length(cause) > 1){
        stop(paste0("Can only predict one cause. Provided are: ", 
                    paste(cause, collapse = ", "), sep = ""))
    }
    if (missing(cause)) {
        cause <- object$theCause
    }
    causes <- object$causes
    if (any(match(as.character(cause), causes, nomatch = 0)==0L))
        stop(paste0("Requested cause ",as.character(cause)," does not match fitted causes which are:\n ",paste0("- ",causes,collapse="\n")))
    ## stopifnot(match(as.character(cause), causes, nomatch = 0) != 
    ## 0)
    if (survtype == "survival") {
        if (object$theCause != cause) 
            stop("Object can be used to predict cause ", object$theCause, 
                 " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
    }
    
    # predict cumulative cause specific hazards
    causeHazard <- predictCox(object$models[[paste("Cause",cause)]],
                              newdata = newdata,
                              times=object$eventTimes,
                              type = c("hazard","cumHazard", if(se){"survival"}else{NULL}), 
                              keep.strata = TRUE,keep.times=TRUE,keep.lastEventTime = TRUE,
                              se = se)
 
    if(keep.strata){
        strata <- causeHazard$strata
    }
    if (survtype == "hazard") {
        cumHazOther <- Reduce("+",lapply(causes[-match(cause, causes)], 
                                         function(c) {
                                             predictCox(object$models[[paste("Cause",c)]],
                                                        newdata = newdata,
                                                        times=object$eventTimes,
                                                        type = "cumHazard",
                                                        keep.strata = FALSE)$cumHazard
                                         }))
        if(se){
          causeCompete <- predictCox(object$models[[paste("Cause",causes[-match(cause, causes)])]],
                                     newdata = newdata,
                                     times=object$eventTimes,
                                     type = "survival",
                                     keep.strata = FALSE,
                                     se = TRUE)
        }
        survProb <- exp(-causeHazard$cumHazard - cumHazOther)
    } else { 
        survProb <- predictCox(object$models[["OverallSurvival"]],
                               type = "survival",
                               times=object$eventTimes,
                               newdata = newdata,
                               keep.strata = FALSE)$survival
    }
    ## system.time(out <- t(apply(survProb * causeHazard$hazard, 1, cumsum)))
    out <- rowCumSum(survProb * causeHazard$hazard)
    ## FIXME: try to get rid of the censored times where nothing happens to F1 earlier
    ## pos <- prodlim::sindex(jump.times = object$eventTimes, eval.times = times)
    pos <- prodlim::sindex(jump.times = causeHazard$times, eval.times = times)
    etimes.max <- max(causeHazard$lastEventTime)
    if(any(times>etimes.max)){ ## Set NA to predictions after the last event 
      pos[times>etimes.max] <- NA
    }
    p <- cbind(0, out)[, pos + 1,drop=FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
                   NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
                   NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    
    if(se){
      outSE <- seCSC(indexTimes = pos,
                   hazardCause1 = causeHazard$hazard, survivalCause1 = causeHazard$survival, survivalCause2 = causeCompete$survival,
                   hazardCause1.se = causeHazard$hazard.se, survivalCause1.se = causeHazard$survival.se, survivalCause2.se = causeCompete$survival.se)
    }
    
    ## Set NA to predictions between the last event of each strata and the last event
    if(!is.null(causeHazard$strata)){
      for(S in names(causeHazard$lastEventTime)){
        test.times <- (times>causeHazard$lastEventTime[S])*(times<=etimes.max)
        if(any(test.times>0)){ 
          newid.S <- causeHazard$strata==S
          p[newid.S,test.times>0] <- NA
          if(se){outSE[newid.S,test.times>0]}
        }
      }
    }
    
    ## export
    if(keep.strata == FALSE && se == FALSE){
      return(p)
    }else{
      out <- list(prediction = p)
      if(keep.strata == TRUE){out$strata <- strata}
      if(se == TRUE){out$prediction.se <- outSE}
      return(out)
    }
}

seCSC <- function(indexTimes,
                  hazardCause1, survivalCause1, survivalCause2,
                  hazardCause1.se, survivalCause1.se, survivalCause2.se){
  
  integrand <- (hazardCause1*survivalCause1*survivalCause2.se)^2 + (hazardCause1*survivalCause1.se*survivalCause2)^2 + (hazardCause1.se*survivalCause1*survivalCause2.se)^2
  out <- c(0, rowCumSum(integrand))[,indexTimes + 1,drop = FALSE] # it is unclear to which value should be set the variance of the absolute risk before the first event. 0???
  return(out)
}
