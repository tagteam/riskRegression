# inverse of the probability of censoring weigths at the subject specific event times
# {{{ method root
subjectWeights <- function(formula,data,method=c("cox","marginal","km","nonpar","aalen","none"),lag=1){
    if (lag!=1 && lag!=0){ stop("lag must be either 0 or 1")}
    method <- tolower(method)
    method <- match.arg(method,c("cox","marginal","km","nonpar","aalen","none"))
    class(method) <- method
    UseMethod("subjectWeights",method)
}
# }}}
# {{{ None: set weights to 1
subjectWeights.none <- function(formula,data,method,lag=1){
    weights <- rep(1,NROW(data))
    out <- list(weights=weights,
                fit=NULL,
                call=match.call(),
                method=method)
    class(out) <- "subjectWeights"
    out
}
subjectWeights.none <- subjectWeights.none
# }}}
# {{{ reverse Kaplan-Meier 
subjectWeights.marginal <- function(formula,data,method,lag=1){
  formula <- update.formula(formula,"~1")
  fit <- prodlim(formula,data=data,reverse=TRUE)
  weights <- predictSurvIndividual(fit,lag=lag)
  out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
  class(out) <- "subjectWeights"
  out
}
subjectWeights.km <- subjectWeights.marginal
# }}}
# {{{ reverse Stone-Beran 
subjectWeights.nonpar <- function(formula,data,method,lag=1){
  fit <- prodlim(formula,data=data,reverse=TRUE,bandwidth="smooth")
  weights <- predictSurvIndividual(fit,lag=lag)
  out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
  class(out) <- "subjectWeights"
  out
}
# }}}
# {{{ reverse Cox via Harrel's package
subjectWeights.cox <- function(formula,data,method,lag=1){
  ## require(rms)
  status.name <- all.vars(formula)[2]
  time.name <- all.vars(formula)[1]
  times <- data[,time.name]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  stopifnot(NROW(na.omit(data))>0)
  fit <- cph(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
  if (lag==1)
    weights <- as.vector(survest(fit,times=times-min(diff(c(0,unique(times))))/2,what='parallel'))
  else # (lag==0)
    weights <- as.vector(survest(fit,times=times,what='parallel'))
  out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
  class(out) <- "subjectWeights"
  out
}

# }}}
# {{{ reverse Aalen method via the timereg package
## subjectWeights.aalen <- function(formula,data,method,lag=1){
  ## require(timereg)
  ## require(rms)
  ## status.name <- all.vars(formula)[2]
  ## time.name <- all.vars(formula)[1]
  ## times <- data[,time.name]
  ## reverse.data <- data
  ## reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  ## fit <- do.call(method,list(formula=formula,data=reverse.data,n.sim=0))
  ## if (lag==1) 
    ## weights <- diag(predictSurvProb(fit,newdata=data,times=pmax(0,times-min(diff(unique(times)))/2)))
  ## else # (lag==0)
  ## weights <- diag(predictSurvProb(fit,newdata=data,times=times))
  ## out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
  ## class(out) <- "subjectWeights"
  ## out
## }
## # }}}
