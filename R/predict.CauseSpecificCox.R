## predict.CauseSpecificCox <- function(object,newdata,times,cause,...){
  ## survtype <- object$survtype
  ## N <- NROW(newdata)
  ## NC <- length(object$model)
  ## eTimes <- object$time
  ## if (missing(cause))
    ## cause <- object$theCause
  ## causes <- object$causes
  ## stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
  ## # predict cumulative cause specific hazards
  ## baseHaz1 <- baselineHazard.coxph(object$models[[paste("Cause",cause)]],times=eTimes)
  ## cumHaz1 <- predict(object$models[[paste("Cause",cause)]],type="risk")%*%baseHaz1
  ## cumHaz1a <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata))
  ## if (length(eTimes)==1)
    ## Haz1 <- cumHaz1
  ## else
    ## Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  ## if (survtype=="hazard"){
    ## cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      ## -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    ## })
    ## lagsurv <- exp(-cumHaz1- do.call("+",cumHazOther))
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## else{
    ## tdiff <- min(diff(eTimes))/2
    ## lagsurv <- predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## pos <- sindex(jump.times=eTimes, eval.times=times)
  ## cbind(0,cuminc1)[,pos+1,drop=FALSE]
## }
