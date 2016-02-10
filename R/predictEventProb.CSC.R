predictEventProb.CSC <- function (object, newdata, times, cause, time, status, Z, ...){
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
    trycumhaz1 <- try(cumHaz1 <- -log(riskRegression:::predictSurvProb.coxph(object$models[[paste("Cause", 
                                                                            cause)]], times = eTimes, newdata = newdata,time,status,Z,cause)), silent = FALSE)
    if (inherits(trycumhaz1, "try-error") == TRUE) 
        stop("Prediction of cause-specific Cox model failed")
    if (length(eTimes) == 1) 
        Haz1 <- cumHaz1
    else Haz1 <- t(apply(cbind(0, cumHaz1), 1, diff))
    if (survtype == "hazard") {
        cumHazOther <- lapply(causes[-match(cause, causes)], 
                              function(c) {
                                  trycumhaz <- try(cumHaz.c <- -log(riskRegression:::predictSurvProb.coxph(object$models[[paste("Cause", 
                                                                                                          c)]], times = eTimes, newdata = newdata, time, status, Z,cause=as.numeric(causes[-match(cause, causes)]))), 
                                                   silent = FALSE)
                                  if (inherits(trycumhaz, "try-error") == TRUE) 
                                      stop("Prediction of cause-specific Cox model failed")
                                  cumHaz.c
                              })
        lagsurv <- exp(-cumHaz1 - Reduce("+", cumHazOther))
        cuminc1 <- t(apply(lagsurv * Haz1, 1, cumsum))
    }
    else {
        tdiff <- min(diff(eTimes))/2
        trylagsurv <- try(lagsurv <- riskRegression:::predictSurvProb.coxph(object$models[["OverallSurvival"]], 
                                                      times = eTimes - tdiff, newdata = newdata, time, status, Z,cause), silent = FALSE)
        if (inherits(trylagsurv, "try-error") == TRUE) 
            stop("Prediction of overall curvival Cox model failed")
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
