oldpredictEventProb.cumIncCox <- function (object, newdata, times, cause, ...) {
  require(survival)
  call <- match.call()
  N <- NROW(newdata)
  # event times for all causes
  csTimes <- object$eventTimes
  overall.haz <- object$survCox
  predSurv <- predictSurvProb(overall.haz,newdata = newdata,times = csTimes)
  laggedSurv <- cbind(rep(1, NROW(newdata)), predSurv)[, 1:length(csTimes),drop=FALSE]
  haz <- object$csHazCox
  cum.bhaz <- basehaz(haz,centered=FALSE)
  cum.bhaz <- c(0, cum.bhaz[, "hazard"])[1 + sindex(jump.times = cum.bhaz[,"time"], eval.times = csTimes)]
  bhaz <- c(cum.bhaz[1], diff(cum.bhaz))
  linPred <- predict(haz, newdata = newdata, type = "lp", se.fit = TRUE)
  cumH <- do.call("cbind", lapply(1:N, function(n){
    cumhaz.n = bhaz * exp(linPred$fit[n, ,drop=FALSE])
  }))
  p <- do.call("rbind", lapply(1:N, function(n) {
    cumsum(laggedSurv[n, ,drop=FALSE] * cumH[, n])
  }))
  p
  pos <- sindex(jump.times=csTimes, eval.times=times)
  cbind(0,p)[,pos+1,drop=FALSE]
}




