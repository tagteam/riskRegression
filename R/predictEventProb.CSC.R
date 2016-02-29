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

##### 
## Questions: why cause=as.numeric(causes[-match(cause, causes)]) indead of c in the lapply in the initial function
##            why test whether predictSurvProb.coxph fails
##            survtype == "hazard" or what

#' @title Predicting hazard or cumulative hazard
#' 
#' @aliases predictHazard predictHazard.coxph predictSurvProb predictSurvProb.coxph
#
#' @param object The fitted coxph model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param cause should the hazard or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' @return a data table containing the predictions for each patient (in rows) and each time (in columns). 
#' 
#' @examples 
#' 
#' data <- prodlim:::SimCompRisk(n)
#' data$time <- round(data$time,1)
#' CSC.fit <- CSC(Hist(time,event)~ X1*X2,data=data, ties = "breslow" )
#' 
#' predCSC <- predictEventProb(CSC.fit, newdata = data, times = unique(sort(data$time)), cause = 1, method = "cpp")
#' 
#' @export
#' 

predictEventProb.CauseSpecificCox <- function (object, newdata, times, cause, method.baseHaz = "cpp"){
  require(pec)
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
                                 newdata = newdata, type = c("hazard","cumHazard"), method = "cpp", col.strata = FALSE)
  
  cumHaz1 <- resPred$cumHazard
  Haz1 <- resPred$hazard
  
  if (survtype == "hazard") {
    cumHazOther <- lapply(causes[-match(cause, causes)], 
                          function(c) {
                            predictHazard.coxph(object$models[[paste("Cause",c)]], times = eTimes, 
                                                newdata = newdata, type = "cumHazard", col.strata = FALSE)
                          })
    lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
    cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
  }
  else {
    tdiff <- min(diff(eTimes))/2
    lagsurv <- predictSurvProb.coxph(object$models[["OverallSurvival"]], times = eTimes - tdiff, newdata = newdata)
    cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
  }
  pos <- prodlim::sindex(jump.times = eTimes, eval.times = times)
  p <- cbind(0, cuminc1)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}
