timevarCoef.riskRegression <- function(object,times,...){
  tvars <- all.vars(object$design$timevar$formula)
  Flevels <- object$factorLevels
  if (length(tvars)==0){
    "none"
  }
  else{
    if (missing(times)) times <- quantile(object$time)
    showTimes <- sindex(eval.times=times,jump.times=object$time)
    showMat <- cbind(times=showTimes,exp(object$timeVarCoef[showTimes,-1,drop=FALSE]))
    rownames(showMat) <- signif(object$timeVarCoef[showTimes,1],2)
    showMat
  }
}
