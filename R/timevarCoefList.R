timevarCoefList <- function(object,times,bytime=TRUE,...){
  if (missing(times)) times <- quantile(object$time)
  showTimes <- sindex(eval.times=times,jump.times=object$time)
  if (bytime==TRUE){
    clist <- lapply(showTimes,function(t){
      coefTable(object$design$timevar$formula,coef=object$timeVarCoef[t,-1,drop=TRUE],se=sqrt(object$timeVarVar[t,-1,drop=TRUE]),levels=object$factorLevels,ref=object$refLevels,trans="exp",test="wald")
    })
  }
  else{

  }    
  names(clist) <- times
  class(clist) <- "timevarCoefList"
  clist
}
