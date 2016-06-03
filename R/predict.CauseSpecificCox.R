#' @title Predicting hazard or cumulative hazard
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
#' dt <- SimCompRisk(1e2)
#' dt$time <- round(dt$time,1)
#' 
#' ttt <- sample(x = unique(sort(dt$time)), size = 10) 
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=dt, method = "breslow" )
#' 
#' predCSC <- predict(CSC.fit, newdata = dt, cause = 2, times = ttt)
#' 
#' #### cph function
#' CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=dt, method = "breslow", fitter = "cph", y = TRUE)
#' 
#' predcph <- predict(CSC.cph, newdata = dt, cause = 2, times = ttt)
#'
#' 
#' 
#' @export
predict.CauseSpecificCox <- function(object,newdata,times,cause,keep.strata = FALSE,...){
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
                              type = c("hazard","cumHazard"), 
                              keep.strata = keep.strata,keep.times=TRUE)
    if(keep.strata){
        strata <- causeHazard$hazard$strata
        causeHazard$hazard[, strata:=NULL]
        causeHazard$cumHazard[, strata := NULL]
    }
    ncol.pred <-  ncol(causeHazard$cumHazard)
    if (survtype == "hazard") {
        cumHazOther <- Reduce("+",lapply(causes[-match(cause, causes)], 
                                         function(c) {
                                             predictCox(object$models[[paste("Cause",c)]],
                                                        newdata = newdata,
                                                        times=object$eventTimes,
                                                        type = "cumHazard",
                                                        keep.strata = FALSE)$cumHazard
                                         }))
        survProb <- exp(-causeHazard$cumHazard - cumHazOther)
    } else { 
        survProb <- predictCox(object$models[["OverallSurvival"]],
                               type = "cumHazard",
                               times=object$eventTimes,
                               newdata = newdata,
                               keep.strata = FALSE)$cumHazard
    }
    ## system.time(out <- t(apply(survProb * causeHazard$hazard, 1, cumsum)))
    out <- rowCumSum(as.matrix(survProb * causeHazard$hazard))
    ## FIXME: try to get rid of the censored times where nothing happens to F1 earlier
    ## pos <- prodlim::sindex(jump.times = object$eventTimes, eval.times = times)
    pos <- prodlim::sindex(jump.times = causeHazard$times, eval.times = times)
    p <- cbind(0, out)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
                   NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
                   NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    ## export
    if(keep.strata == TRUE){
        return(list(p, strata = strata))
    }else{
        p
    }
}
