#' @title Compute the baseline hazard
#
#' @param object The fitted coxph model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated probabilities
#' @param cause The event of interest
#' @param tYpe should the survival or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' We need to decide what to do after the last event
#' Strange behavior of predict when strata
#' What to return if stratified Cox model with no covariate
#' Can use copy to have direct display of the result
#' 
#' @return a data table containing the predictions for each patient (in rows) and each time (in columns) and the strata (if any).
#' 
#' @examples 
#' 
#' library(data.table)
#' library(survival)
#' library(prodlim)
#' library(pec)    # needed to define the  predictSurvProb method
#' 
#' set.seed(10)
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#'
#' res1 <- predictSurvProb(fit, newdata = d, times = 10)
#' res2 <- predictSurvProb(fit, newdata = d, times = d$time)
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- predictSurvProb(fitS, newdata = d, times = 10)
#' res2S <- predictSurvProb(fitS, newdata = d, times = d$time)
#' 
#' @export
#' 
predictSurvProb.coxph <- function (object, newdata, times, cause = 1, type = "survival",
                                   method.baseHaz = "dt") {
  
  match.arg(type, choices = c("hazard", "cumHazard","survival"), several.ok = FALSE)
  
  times <- sort(times)
  Lambda0 <- baseHaz(object, cause = cause, method = method.baseHaz, centered = TRUE, lasttime =  times[length(times)],
                     addFirst = TRUE, addLast = TRUE) 
  
  strataspecials <- attr(object$terms,"specials")$strata 

  if(is.null(strataspecials)){
    Xb <- predict(object, newdata, type = "lp")
    
    ## compute survival
    time.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times)
    
    resPred <- switch(type,
                      "hazard" = data.table( exp(Xb) %o% Lambda0$hazard[time.index] ),
                      "cumHazard" = data.table( exp(Xb) %o% Lambda0$cumHazard[time.index] ),
                      "survival" = data.table( exp(- exp(Xb) %o% Lambda0$cumHazard[time.index]) )
                      
    )
    
    setnames(resPred, paste0("t",times))
    
  }else{
    
    if("XXstrata" %in% names(newdata)){
      stop("predictSurvProb2.coxph: \'newdata\' must not contains a column named \"XXstrata\" \n")
    }
      
    require(data.table)
    newdata <- as.data.table(newdata)
    
    ## linear predictor
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    Terms <- attr(object$terms, "term.labels")
    BaseVar <- Terms[ Terms %in% stratalabels == FALSE]
    if(length(BaseVar)>0){ # predict crashes if no additional variable
      Xb <-  rowSums(predict(object, newdata, type = "terms")) 
      # predict(object, newdata, type = "lp") has an output that I do not understand and that does not match "terms" in presence of strata
    }else{
      Xb <- rep(0, nrow(newdata))
    }
    
    ## define the strata
    stratalabels <- attr(object$terms,"term.labels")[strataspecials - 1]
    names.newdata<- names(newdata)
    strataVar <- names.newdata[which(paste0("strata(", names.newdata,")") %in% stratalabels)]
    nStrata <- length(strataVar)
    sapply(1:nStrata, function(x){
      newdata[, strataVar[x] := paste0(strataVar[x],"=",.SD[[1]]),.SDcols = strataVar[x], with = FALSE]
    })
    newdata[, XXstrata := interaction(.SD, drop = TRUE, sep = ", ") ,.SDcols = strataVar]
    
    ## compute survival
    resPred <- NULL
    levelsStrata <- levels(newdata$XXstrata)
    for(iterS in levelsStrata){
      indexHaz <- Lambda0[,.I[strata == iterS]]
      indexNew <- newdata[,.I[XXstrata == iterS]]
      time.index <- prodlim::sindex(jump.times=Lambda0[indexHaz,time],eval.times=times)  
      resPred_tempo <- switch(type,
                              "hazard" = exp(Xb[indexNew]) %o% Lambda0[indexHaz[time.index],hazard],
                              "cumHazard" = exp(Xb[indexNew]) %o% Lambda0[indexHaz[time.index],cumHazard],
                              "survival" = exp(- exp(Xb[indexNew]) %o% Lambda0[indexHaz[time.index],cumHazard] )
      )
      
      resPred <-  rbindlist(list(resPred, 
                                 data.table(resPred_tempo, strata = iterS, index = indexNew))
                            )
    }
    setkey(resPred,index)
    resPred[,index := NULL]
    setnames(resPred, c(paste0("t",times),"strata"))
    
  }
  
  # export
  return(resPred)
}



#
# Bazeline Hazard (Cum.)
#

# h_theta_breslow_strata <- function(object,time, status, Z,cause){
#   
#   strataspecials = attr(object$terms,"specials")$strata
#   if (is.null(strataspecials)) {
#     eventtimes <- unique(sort(time[status!=0]))
#     h_theta_event <- numeric(length=length(unique(eventtimes)))
#     coef <- object$coef[sort(names(object$coef))]
#     
#     for (i in order(eventtimes)){
#       Ri <- (time>=eventtimes[i])
#       ExpCov <- matrix(NA, nrow=sum(Ri),ncol=1)
#       for (j in 1:sum(Ri)) { 
#         ExpCov[j] <- exp(coef%*%t(Z[Ri,])[,j])
#       }
#       h_theta_event[i] <- sum((time==eventtimes[i]) & (status==cause))/sum(ExpCov)  
#     }
#     
#     h_theta <- h_theta_event[match(unique(sort(time)),eventtimes)]
#     h_theta[is.na(h_theta)] <- 0
#     H_theta <- cumsum(h_theta)
#     new_basehaz <- data.frame(hazard = H_theta, time = sort(unique(time)))
#     return(new_basehaz)
#   } 
#   else {
#     new_basehaz <- NULL
#     stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
#     strata <- object$xlevels[[stratalabels]] 
#     mod <- model.frame(object) 
#     Terms <- attr(object$terms, "term.labels")
#     BaseVar <- Terms[ Terms != stratalabels ]
#     
#     for (s in 1:length(strata)){
#       Si <- mod[stratalabels][,1]==levels(mod[stratalabels][,1])[s]
#       strataeventtimes <- unique(sort(time[status!=0 & Si]))
#       stratatime <- time[Si]
#       h_theta_event <- numeric(length=length(unique(strataeventtimes)))
#       for (i in order(strataeventtimes)){
#         Ri <- (stratatime>=strataeventtimes[i])
#         ExpCov <- matrix(NA, nrow=sum(Ri),ncol=1)
#         for (j in 1:sum(Ri)) { 
#           if(length(BaseVar)!=1){
#             ExpCov[j] <- exp(object$coef%*%t((Z[Si,BaseVar][Ri,][j,])))}
#           else {
#             ExpCov[j] <- exp(object$coef%*%t((Z[Si,BaseVar][Ri][j])))
#           }
#         }
#         h_theta_event[i] <- sum((stratatime==strataeventtimes[i]) & (status[Si]==cause))/sum(ExpCov)
#       }
#       h_theta <- h_theta_event[match(unique(sort(stratatime)),strataeventtimes)]
#       h_theta[is.na(h_theta)] <- 0
#       H_theta <- cumsum(h_theta)
#       new_basehaz <- rbind(new_basehaz, data.frame(hazard = H_theta,
#                                                    time = sort(unique(stratatime)), strata=levels(mod[stratalabels][,1])[s]))
#     }
#     return(new_basehaz)
#   }
# }

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
