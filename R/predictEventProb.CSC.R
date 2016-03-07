#' @export
predictEventProb <- function (x, ...) {
  UseMethod("predictEventProb", x)
}

# predictEventProb.CSC <- function (object, newdata, times, cause, time, status, Z, ...){
#   require(pec)
#     survtype <- object$survtype
#     N <- NROW(newdata)
#     NC <- length(object$model)
#     if (length(cause) > 1) 
#         stop(paste0("Can only predict one cause. Provided are: ", 
#                     paste(cause, collapse = ", "), sep = ""))
#     eTimes <- object$eventTimes
#     if (missing(cause)) 
#         cause <- object$theCause
#     causes <- object$causes
#     stopifnot(match(as.character(cause), causes, nomatch = 0) != 
#                   0)
#     if (survtype == "survival") {
#         if (object$theCause != cause) 
#             stop("Object can be used to predict cause ", object$theCause, 
#                  " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
#     }
#     trycumhaz1 <- try(cumHaz1 <- -log(riskRegression:::predictSurvProb.coxph(object$models[[paste("Cause", 
#                                                                             cause)]], times = eTimes, newdata = newdata,time,status,Z,cause)), silent = FALSE)
#     if (inherits(trycumhaz1, "try-error") == TRUE) 
#         stop("Prediction of cause-specific Cox model failed")
#     if (length(eTimes) == 1) 
#         Haz1 <- cumHaz1
#     else Haz1 <- t(apply(cbind(0, cumHaz1), 1, diff))
#     if (survtype == "hazard") {
#         cumHazOther <- lapply(causes[-match(cause, causes)], 
#                               function(c) {
#                                   trycumhaz <- try(cumHaz.c <- -log(riskRegression:::predictSurvProb.coxph(object$models[[paste("Cause", 
#                                                                                                           c)]], times = eTimes, newdata = newdata, time, status, Z,cause=as.numeric(causes[-match(cause, causes)]))), 
#                                                    silent = FALSE)
#                                   if (inherits(trycumhaz, "try-error") == TRUE) 
#                                       stop("Prediction of cause-specific Cox model failed")
#                                   cumHaz.c
#                               })
#         lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
#         cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
#     }
#     else {
#         tdiff <- min(diff(eTimes))/2
#         trylagsurv <- try(lagsurv <- riskRegression:::predictSurvProb.coxph(object$models[["OverallSurvival"]], 
#                                                       times = eTimes - tdiff, newdata = newdata, time, status, Z,cause), silent = FALSE)
#         if (inherits(trylagsurv, "try-error") == TRUE) 
#             stop("Prediction of overall curvival Cox model failed")
#         cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
#     }
#     pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
#     p <- cbind(0, cuminc1)[, pos + 1, drop = FALSE]
#     if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
#         stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
#                    NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
#                    NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
#     p
# }

#' @title Predicting hazard or cumulative hazard
#' 
#' @aliases predictEventProb predictEventProb.CauseSpecificCox
#
#' @param object The fitted coxph model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param cause should the hazard or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' 
#' @details 
#' Not suited for time varying cox models.
#' 
#' @return A data table containing the predictions for each patient (in rows) and each time (in columns). 
#' 
#' @examples 
#' 
#' data <- SimCompRisk(1e3)
#' data$time <- round(data$time,1)
#' CSC.fit <- CSC(Hist(time,event)~ X1*X2,data=data, ties = "breslow" )
#' 
#' predCSC <- predictEventProb(CSC.fit, newdata = data, cause = 2, times = unique(sort(data$time)))
#' 
#' @export
predictEventProb.CauseSpecificCox <- function (object, newdata, times, cause,
                                               col.strata = FALSE, method.baseHaz = "cpp"){
  survtype <- object$survtype
  N <- NROW(newdata)
  NC <- length(object$model)
  if (length(cause) > 1) 
    stop(paste0("Can only predict one cause. Provided are: ", 
                paste(cause, collapse = ", "), sep = ""))
  eTimes <- object$eventTimes
  if (missing(cause)) 
    cause <- object$theCause
  causes <- object$causes
  stopifnot(match(as.character(cause), causes, nomatch = 0) != 
              0)
  if (survtype == "survival") {
    if (object$theCause != cause) 
      stop("Object can be used to predict cause ", object$theCause, 
           " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
  }
  
  resPred <- predictHazard.coxph(object$models[[paste("Cause",cause)]], times = eTimes, 
                                 newdata = newdata, type = c("hazard","cumHazard"), method = method.baseHaz, col.strata = col.strata)
  
  
  range(t(apply(resPred[[2]], 1, diff)) - resPred[[1]][,-1, with = FALSE], na.rm = TRUE)
 
  if(col.strata){
    strata <- resPred$hazard$strata
    resPred$hazard[, strata := NULL]
    resPred$cumHazard[, strata := NULL]
  }
  ncol.pred <-  ncol(resPred$cumHazard)
  
  if (survtype == "hazard") {
    cumHazOther <- lapply(causes[-match(cause, causes)], 
                          function(c) {
                            predictHazard.coxph(object$models[[paste("Cause",c)]], times = eTimes, 
                                                newdata = newdata, type = "cumHazard", method = method.baseHaz, col.strata = FALSE)
                          })

    resPred$cumHazard[, 1:ncol.pred := resPred$hazard * exp(-.SD - Reduce("+", cumHazOther)) , .SD = 1:ncol.pred]
    
    resPred$cumHazard[, 1:ncol.pred := unlist(apply(apply(.SD, 1, cumsum), 1 , list), recursive = FALSE),
                      .SD = 1:ncol.pred] 
    
    resPred$cumHazard[, tFirst := 0]
    
    data.table::setcolorder(resPred$cumHazard, c("tFirst", setdiff(names(resPred$cumHazard), "tFirst")))
    
  } else { ### NOT TESTED
    
    tdiff <- min(diff(eTimes))/2
    
    cuminc1 <- predictSurvProb.coxph(object$models[["OverallSurvival"]], times = eTimes - tdiff, 
                                   newdata = newdata, method = method.baseHaz, col.strata = FALSE)
    
    cuminc1[, 1:ncol.pred := unlist(apply(apply(.SD* Haz1, 1, cumsum), 1 , list), recursive = FALSE),
                  .SD = 1:ncol.pred] 
  }

  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- resPred$cumHazard[, pos + 1, with = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  
  ## export
  if(col.strata == TRUE){
    return(data.table::data.table(p, strata = strata))
  }else{
    p
  }
}
