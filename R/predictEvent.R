#' @title Predicting hazard or cumulative hazard
#' 
#' @aliases predictEvent
#
#' @param object The fitted cause specific Cox model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param cause should the hazard or the cumulative hazard be returned
#' @param keep.strata If \code{TRUE} return an additional column containing the strata level. 
#' 
#' @details 
#' Note: for Cox regression models with time varying covariates it does not make sense to use this function, 
#' because the predicted risk has to be a measurable function of the data available at the time origin. 
#' 
#' @return A data table containing the predictions for each patient (in rows) and each time (in columns). 
#' 
#' @examples 
#' 
#' data <- SimCompRisk(1e3)
#' data$time <- round(data$time,1)
#' 
#' seq_times <- sample(x = unique(sort(data$time)), size = 100) 
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1*X2,data=data, ties = "breslow" )
#' 
#' predCSC <- predictEvent(CSC.fit, newdata = data, cause = 2, times = seq_times)
#' 
#' #### cph function
#' CSC.fit <- CSC(Hist(time,event)~ X1*X2,data=data, ties = "breslow", fitter = "cph", y = TRUE)
#' 
#' predCSC <- predictEvent(CSC.fit, newdata = data, cause = 2, times = seq_times)
#' @export
predictEvent <- function (object, newdata, times, cause,
                          keep.strata = FALSE){
  survtype <- object$survtype
  if (length(cause) > 1){
    stop(paste0("Can only predict one cause. Provided are: ", 
                paste(cause, collapse = ", "), sep = ""))
  }
  eTimes <- object$eventTimes
  if (missing(cause)) {
    cause <- object$theCause
  }
  causes <- object$causes
  stopifnot(match(as.character(cause), causes, nomatch = 0) != 
              0)
  if (survtype == "survival") {
    if (object$theCause != cause) 
      stop("Object can be used to predict cause ", object$theCause, 
           " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
  }
  
  resPred <- predictSurv(object$models[[paste("Cause",cause)]], times = eTimes, 
                         newdata = newdata, type = c("hazard","cumHazard"), 
                         format = if(keep.strata){"data.table"}else{"matrix"}, keep.strata = keep.strata)
  
  if(keep.strata){
    strata <- resPred$hazard$strata
    resPred$hazard[, strata := NULL]
    resPred$cumHazard[, strata := NULL]
  }
  ncol.pred <-  ncol(resPred$cumHazard)
  
  if (survtype == "hazard") {
    cumHazOther <- lapply(causes[-match(cause, causes)], 
                          function(c) {
                            predictSurv(object$models[[paste("Cause",c)]], times = eTimes, 
                                        newdata = newdata, type = "cumHazard", 
                                        format = "matrix", keep.strata = FALSE)
                          })
  
    resPred <- t(apply(resPred$hazard * exp(-resPred$cumHazard - Reduce("+", cumHazOther)), 1, cumsum))
    
  } else { ### NOT TESTED
    
    tdiff <- min(diff(eTimes))/2
    
    resPred <- predictSurv(object$models[["OverallSurvival"]], times = eTimes - tdiff, type = "cumHazard", 
                           newdata = newdata, format = "matrix", keep.strata = FALSE)
    
    resPred <- t(apply(resPred * Haz1, 1, cumsum))
  }
  
  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- cbind(0, resPred)[, pos + 1]
  
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  
  ## export
  if(keep.strata == TRUE){
    return(data.table::data.table(p, strata = strata))
  }else{
    data.table::data.table(p)
    # p
  }
}
