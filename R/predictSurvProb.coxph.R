#
# Bazeline Hazard (Cum.)
#

h_theta_breslow_strata <- function(object,time, status, Z,cause){
  
  strataspecials = attr(object$terms,"specials")$strata
  if (is.null(strataspecials)) {
    eventtimes <- unique(sort(time[status!=0]))
    h_theta_event <- numeric(length=length(unique(eventtimes)))
    coef <- object$coef[sort(names(object$coef))]
    
    for (i in order(eventtimes)){
      Ri <- (time>=eventtimes[i])
      ExpCov <- matrix(NA, nrow=sum(Ri),ncol=1)
      for (j in 1:sum(Ri)) { 
        ExpCov[j] <- exp(coef%*%t(Z[Ri,])[,j])
      }
      h_theta_event[i] <- sum((time==eventtimes[i]) & (status==cause))/sum(ExpCov)  
    }
    
    h_theta <- h_theta_event[match(unique(sort(time)),eventtimes)]
    h_theta[is.na(h_theta)] <- 0
    H_theta <- cumsum(h_theta)
    new_basehaz <- data.frame(hazard = H_theta, time = sort(unique(time)))
    return(new_basehaz)
  } 
  else {
    new_basehaz <- NULL
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    strata <- object$xlevels[[stratalabels]] 
    mod <- model.frame(object) 
    Terms <- attr(object$terms, "term.labels")
    BaseVar <- Terms[ Terms != stratalabels ]
    
    for (s in 1:length(strata)){
      Si <- mod[stratalabels][,1]==levels(mod[stratalabels][,1])[s]
      strataeventtimes <- unique(sort(time[status!=0 & Si]))
      stratatime <- time[Si]
      h_theta_event <- numeric(length=length(unique(strataeventtimes)))
      for (i in order(strataeventtimes)){
        Ri <- (stratatime>=strataeventtimes[i])
        ExpCov <- matrix(NA, nrow=sum(Ri),ncol=1)
        for (j in 1:sum(Ri)) { 
          if(length(BaseVar)!=1){
            ExpCov[j] <- exp(object$coef%*%t((Z[Si,BaseVar][Ri,][j,])))}
          else {
            ExpCov[j] <- exp(object$coef%*%t((Z[Si,BaseVar][Ri][j])))
          }
        }
        h_theta_event[i] <- sum((stratatime==strataeventtimes[i]) & (status[Si]==cause))/sum(ExpCov)
      }
      h_theta <- h_theta_event[match(unique(sort(stratatime)),strataeventtimes)]
      h_theta[is.na(h_theta)] <- 0
      H_theta <- cumsum(h_theta)
      new_basehaz <- rbind(new_basehaz, data.frame(hazard = H_theta,
                                                   time = sort(unique(stratatime)), strata=levels(mod[stratalabels][,1])[s]))
    }
    return(new_basehaz)
  }
}

#
# Own PredictSurvProb
#

# predictSurvProb.coxph<- function (object, times, newdata, time, status, Z, cause) {
#   NewBasehaz <- h_theta_breslow_strata(object, time,status, Z,cause)
#   p <- matrix(NA, nrow=nrow(newdata), ncol=length(times))
#   strataspecials = attr(object$terms,"specials")$strata
#   
#   if (is.null(strataspecials)){
#     time.index = sort(sindex(jump.times=NewBasehaz$time,eval.times=times))
#     for (i in 1:nrow(newdata)) {
#       for (j in 1:length(times)) { 
#         t <-  sort(times)[j]
#         p[i,j] <- exp(-c(0,NewBasehaz[,1])[1+time.index[j]])^exp(sum(newdata[i,names(Z)]*object$coef))
#       }
#     }
#     return(p)
#   }
#   else{
#     stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
#     Terms <- attr(object$terms, "term.labels")
#     BaseVar <- Terms[ Terms != stratalabels ]
#     
#     for (i in 1:nrow(newdata)) {
#       Base <- NewBasehaz[NewBasehaz$strata==eval(parse(text = stratalabels), newdata)[i],]
#       newdata$strata <- eval(parse(text = stratalabels), newdata)
#       time.index = sort(sindex(jump.times=Base$time,eval.times=times))
#       
#       for (j in 1:length(times)) { 
#         t <-  sort(times)[j]   
#         Si <- (newdata$strata == eval(parse(text = stratalabels), newdata)[i])
#         p[i,j] <- exp(-c(0,Base[,1])[1+time.index[j]])^exp(sum(newdata[i,BaseVar]*object$coef))
#       }
#     }
#     
#     return(p)
#   }
# }

predictSurvProb.coxph <- function (object, times, newdata, cause = 1, method.baseHaz = "dt") {
  
  times <- sort(times)
  Lambda0 <- baseHaz(object, cause = cause, method = method.baseHaz, center = FALSE, lasttime = times[length(times)])
  Xb <- predict(object, newdata, type = "lp") # seems to crash if only strata variables
  
  strataspecials <- attr(object$terms,"specials")$strata 
  
  if(is.null(strataspecials)){
    time.index <- prodlim::sindex(jump.times = c(-1, Lambda0$time), eval.times = times) # -1 if times before first event
    resSurv <- exp(exp(Xb) %o% c(0,-Lambda0$cumHazard)[time.index]) 
    resSurv <- data.table( resSurv )
    setnames(resSurv, paste0("t",times))
    
  }else{
    
    if("XXstrata" %in% names(newdata)){
      stop("predictSurvProb2.coxph: \'newdata\' must not contains a column named \"XXstrata\" \n")
    }
      
    require(data.table)
    newdata <- as.data.table(newdata)
    
    # define the strata
    stratalabels <- attr(object$terms,"term.labels")[strataspecials - 1]
    names.newdata<- names(newdata)
    strataVar <- names.newdata[which(paste0("strata(", names.newdata,")") %in% stratalabels)]
    nStrata <- length(strataVar)
    sapply(1:nStrata, function(x){
      newdata[, strataVar[x] := paste0(strataVar[x],"=",.SD[[1]]),.SDcols = strataVar[x], with = FALSE]
    })
    newdata[, XXstrata := interaction(.SD, drop = TRUE, sep = ", ") ,.SDcols = strataVar]
    
     # compute the survival
    resSurv <- NULL
    levelsStrata <- levels(newdata$XXstrata)
    for(iterS in levelsStrata){
      indexTempo <- Lambda0[,.I[strata == iterS]]
      time.index <- prodlim::sindex(jump.times=c(-1,Lambda0[indexTempo,time]),eval.times=times)  
      
      resSurv_tempo <- exp(exp(Xb[newdata[, I(XXstrata == iterS)]]) %o% c(0,-Lambda0$cumHazard[indexTempo])[time.index] )
      resSurv <-  rbindlist(list(resSurv, 
                                 data.table(resSurv_tempo, iterS)),
                            use.names = TRUE, fill = TRUE)
    }
    setnames(resSurv, c(paste0("t",times),"strata"))
    
  }
  
  # export
  return(resSurv)
}