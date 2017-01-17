### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02) 
## Version: 
## last-updated: Oct 31 2016 (07:38) 
##           By: Thomas Alexander Gerds
##     Update #: 37
#----------------------------------------------------------------------
## 
### Commentary:
#
# methods for extracting predictions from various objects
# 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# --------------------------------------------------------------------
#' Extrating predicting risks from regression models 
#' 
#' Extract event probabilities from fitted regression models and machine learning objects. 
#' 
#' The function predictRisk is a generic function, meaning that it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' See \code{\link{predictRisk}}.
#' 
#' @aliases predictRisk predictRisk.CauseSpecificCox
#' predictRisk.riskRegression predictRisk.FGR
#' predictRisk.prodlim predictRisk.rfsrc predictRisk.aalen
#' predictRisk.riskRegression predictRisk.cox.aalen
#' predictRisk.coxph predictRisk.cph predictRisk.default
#' predictRisk.rfsrc predictRisk.matrix predictRisk.pecCtree
#' predictRisk.pecCforest predictRisk.prodlim predictRisk.psm
#' predictRisk.selectCox predictRisk.survfit predictRisk.randomForest
#' predictRisk.lrm predictRisk.default predictRisk.glm
#' predictRisk.rpart
#' @usage
#' \method{predictRisk}{glm}(object,newdata,...)
#' \method{predictRisk}{cox.aalen}(object,newdata,times,...)
#' \method{predictRisk}{cph}(object,newdata,times,se=FALSE,iid=FALSE,...)
#' \method{predictRisk}{coxph}(object,newdata,times,se=FALSE,iid=FALSE,...)
#' \method{predictRisk}{matrix}(object,newdata,times,cause,...)
#' \method{predictRisk}{selectCox}(object,newdata,times,...)
#' \method{predictRisk}{psm}(object,newdata,times,...)
#' \method{predictRisk}{survfit}(object,newdata,times,...)
#' \method{predictRisk}{riskRegression}(object,newdata,times,cause,...)
#' \method{predictRisk}{prodlim}(object,newdata,times,cause,...)
#' \method{predictRisk}{rfsrc}(object,newdata,times,cause,...)
#' \method{predictRisk}{FGR}(object,newdata,times,cause,...)
#' \method{predictRisk}{CauseSpecificCox}(object,newdata,times,cause,se=FALSE,iid=FALSE,...)
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output.
#' @param iid Logical. If \code{TRUE} add the influence function corresponding ot the output.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return For binary outcome a vector with predicted risks. For survival outcome with and without
#' competing risks
#' a matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry is a probability and in
#' rows the values should be increasing.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso See \code{\link{predictRisk}}.
#' @keywords survival
#' @examples
#' ## binary outcome
#' library(rms)
#' set.seed(7)
#' x <- abs(rnorm(20))
#' d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
#' nd <- data.frame(y=rbinom(8,1,x/max(x)),x=abs(rnorm(8)),z=rnorm(8))
#' fit <- lrm(y~x+z,d)
#' predictRisk(fit,newdata=nd)
#'
#' ## survival outcome
#' # generate survival data
##' library(prodlim)
##' set.seed(100)
##' d <- sampleData(100,outcome="survival")
##' # then fit a Cox model
##' library(rms)
##' cphmodel <- cph(Surv(time,event)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
##' # or via survival
##' library(survival)
##' coxphmodel <- coxph(Surv(time,event)~X1+X2,data=d,x=TRUE,y=TRUE)
##' 
##' # Extract predicted survival probabilities 
##' # at selected time-points:
##' ttt <- quantile(d$time)
##' # for selected predictor values:
##' ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
##' # as follows
##' predictRisk(cphmodel,newdata=ndat,times=ttt)
##' predictRisk(coxphmodel,newdata=ndat,times=ttt)
##' 
##' # stratified cox model
##' sfit <- coxph(Surv(time,event)~strata(X1)+X2,data=d,x=TRUE,y=TRUE)
##' predictRisk(sfit,newdata=d[1:3,],times=c(1,3,5,10))
##' 
##' ## simulate learning and validation data
##' learndat <- sampleData(100,outcome="survival")
##' valdat <- sampleData(100,outcome="survival")
##' ## use the learning data to fit a Cox model
##' library(survival)
##' fitCox <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
##' ## suppose we want to predict the survival probabilities for all patients
##' ## in the validation data at the following time points:
##' ## 0, 12, 24, 36, 48, 60
##' psurv <- predictRisk(fitCox,newdata=valdat,times=seq(0,60,12))
##' ## This is a matrix with survival probabilities
##' ## one column for each of the 5 time points
##' ## one row for each validation set individual
##' 
##' # Do the same for a randomSurvivalForest model
##' library(randomForestSRC)
##' rsfmodel <- rfsrc(Surv(time,event)~X1+X2,data=learndat)
##' prsfsurv=predictRisk(rsfmodel,newdata=valdat,times=seq(0,60,12))
##' # plot(psurv,prsfsurv)
##' 
##' ## Cox with ridge option
##' f1 <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
##' f2 <- coxph(Surv(time,event)~ridge(X1)+ridge(X2),data=learndat,x=TRUE,y=TRUE)
##' \dontrun{
##' plot(predictRisk(f1,newdata=valdat,times=10),
##'      riskRegression:::predictRisk.coxph(f2,newdata=valdat,times=10),
##'      xlim=c(0,1),
##'      ylim=c(0,1),
##'      xlab="Unpenalized predicted survival chance at 10",
##'      ylab="Ridge predicted survival chance at 10")
##'}
##' 
#' ## competing risks
#' 
#' library(survival)
#' library(riskRegression)
#' library(prodlim)
#' train <- SimCompRisk(100)
#' test <- SimCompRisk(10)
#' cox.fit  <- CSC(Hist(time,cause)~X1+X2,data=train)
#' predictRisk(cox.fit,newdata=test,times=seq(1:10),cause=1)
#'
#' ## with strata
#' cox.fit2  <- CSC(list(Hist(time,cause)~strata(X1)+X2,Hist(time,cause)~X1+X2),data=train)
#' predictRisk(cox.fit2,newdata=test,times=seq(1:10),cause=1)
#' 
#' @export 
predictRisk <- function(object,newdata,...){
  UseMethod("predictRisk",object)
}

##' @export 
predictRisk.default <- function(object,newdata,times,cause,...){
  stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
}

##' @export
predictRisk.double <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.integer <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.numeric <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

##' @export
predictRisk.glm <- function(object,newdata,...){
    if (object$family$family=="binomial")
        return(as.numeric(stats::predict(object,newdata=newdata,type="response")))
    else{ stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
      }
}

##' @export
predictRisk.formula <- function(object,newdata,...){
    ff <- update.formula(object,"NULL~.")
    if (length(all.vars(ff))==1){
        p <- stats::model.frame(ff,newdata)[[1]]
        p
    } else{
        fit <- glm(object,data=newdata,family="binomial")
        predictRisk(fit,newdata=newdata,...)
    }
}

##' @export
predictRisk.BinaryTree <- function(object,newdata,...){
    treeresponse <- party::treeresponse
    sapply(treeresponse(object,newdata=newdata),function(x)x[1])
}

##' @export
predictRisk.lrm <- function(object,newdata,...){
  as.numeric(stats::predict(object,newdata=newdata,type="fitted"))
}

##' @export
predictRisk.rpart <- function(object,newdata,...){
  p <- as.numeric(stats::predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  p
}

##' @export
predictRisk.randomForest <- function(object,newdata,...){
  as.numeric(stats::predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
}

##' @export
predictRisk.matrix <- function(object,newdata,times,cause,...){
    if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
        stop(paste("Prediction matrix has wrong dimensions: ",
                   NROW(object),
                   " rows and ",
                   NCOL(object),
                   " columns.\n But requested are predicted probabilities for ",
                   NROW(newdata),
                   " subjects (rows) in newdata and ",
                   length(times),
                   " time points (columns)",
                   sep=""))
    }
    object
}


predictRisk.aalen <- function(object,newdata,times,...){
    ## require(timereg)
    stop("FIXME")
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    covanames <- names(time.coef)[-(1:2)]
    notfound <- match(covanames,names(newdata),nomatch=0)==0
    if (any(notfound))
        stop("\nThe following predictor variables:\n\n",
             paste(covanames[notfound],collapse=","),
             "\n\nwere not found in newdata, which only provides the following variables:\n\n",
             paste(names(newdata),collapse=","),
             "\n\n")
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    ## hazard <- .C("survest_cox_aalen",
                 ## timehazard=double(ntime*nobs),
                 ## as.double(unlist(time.coef[,-1])),
                 ## as.double(unlist(time.vars)),
                 ## as.integer(ntimevars+1),
                 ## as.integer(nobs),
                 ## as.integer(ntime),PACKAGE="pec")$timehazard
    hazard <- matrix(hazard,ncol=ntime,nrow=nobs,dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-hazard),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(jump.times=objecttime,eval.times=times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

predictRisk.cox.aalen <- function(object,newdata,times,...){
    #  require(timereg)
    ##  The time-constant effects first
    stop("FIXME")
    const <- c(object$gamma)
    names(const) <- substr(dimnames(object$gamma)[[1]],6,nchar(dimnames(object$gamma)[[1]])-1)
    constant.part <- t(newdata[,names(const)])*const
    constant.part <- exp(colSums(constant.part))
    ##  Then extract the time-varying effects
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    ## time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),PACKAGE="pec")$timehazard
    time.part <- matrix(time.part,ncol=ntime,nrow=nobs)
    ## dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-time.part*constant.part),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(objecttime,times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

    
##' @export
predictRisk.coxph <- function(object,newdata,times,se=FALSE,iid=FALSE,...){
    res <- predictCox(object=object,newdata=newdata,times=times, se = se, iid = iid,keep.times=FALSE,keep.lastEventTime=FALSE,type="survival")

    p <- 1-res$survival
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    if(se){attr(p,"se") <- res$survival.se}
    if(iid){attr(p,"iid") <- res$survival.iid}
    return(p)
}
## predictRisk.coxph <- function(object,newdata,times,...){
## baselineHazard.coxph(object,times)
## require(survival)
## new feature of the survival package requires that the
## original data are included
## survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
## survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
## browser(skipCalls=1)
## survfit.object <- survival::survfit(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
## if (is.null(attr(object$terms,"specials")$strata)){
## ## case without strata 
## inflated.pred <- summary(survfit.object,times=times)$surv
## p <- t(inflated.pred)        
## } else{
## ## case with strata 
## inflated.pred <- summary(survfit.object,times=times)
## plist <- split(inflated.pred$surv,inflated.pred$strata)
## p <- do.call("rbind",lapply(plist,function(x){
## beyond <- length(times)-length(x)
## c(x,rep(NA,beyond))
## }))
## ## p <- matrix(inflated.pred,ncol=length(times),byrow=TRUE)
## }
## if ((miss.time <- (length(times) - NCOL(p)))>0)
## p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
## stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
## 1-p
## }

##' @export
predictRisk.coxph.penal <- function(object,newdata,times,...){
  frailhistory <- object$history$'frailty(cluster)'$history
  frailVar <- frailhistory[NROW(frailhistory),1]
  linearPred <- predict(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  basehaz <- basehaz(object)
  bhTimes <- basehaz[,2]
  bhValues <- basehaz[,1]
  survPred <- do.call("rbind",lapply(1:NROW(newdata),function(i){
    (1+frailVar*bhValues*exp(linearPred[i]))^{-1/frailVar}
  }))
  where <- prodlim::sindex(jump.times=bhTimes,eval.times=times)
  p <- cbind(1,survPred)[,where+1]
  if ((miss.time <- (length(times) - NCOL(p)))>0)
    p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  1-p
}


##' @export 
predictRisk.cph <- function(object,newdata,times,se=FALSE,iid=FALSE,...){
    ## if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
    ## p <- rms::survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
    ## if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
    res <- predictCox(object=object,newdata=newdata,times=times,se = se, iid = iid,keep.times=FALSE,keep.lastEventTime=FALSE,type="survival")
    
    p <- 1-res$survival
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    
    if(se){attr(p,"se") <- res$survival.se}
    if(iid){attr(p,"iid") <- res$survival.iid}
    return(p)
}

##' @export
predictRisk.selectCox <- function(object,newdata,times,...){
    predictRisk(object[[1]],newdata=newdata,times=times,...)
}


##' @export 
predictRisk.prodlim <- function(object,newdata,times,cause,...){
    ## require(prodlim)
    if (object$model=="competing.risks" && missing(cause))
        stop(paste0("Cause is missing. Should be one of the following values: ",paste(attr(object$model.response,"states"),collapse=", ")))
    p <- predict(object=object,cause=cause,type="cuminc",newdata=newdata,times=times,mode="matrix",level.chaos=1)
    ## if the model has no covariates
    ## then all cases get the same prediction
    ## in this exceptional case we return a vector
    if (NROW(p)==1 && NROW(newdata)>=1)
        p <- as.vector(p)
    ## p[is.na(p)] <- 0
    if (is.null(dim(p))){
        if (length(p)!=length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    else{
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    p
}

predict.survfit <- function(object,newdata,times,bytimes=TRUE,type="cuminc",fill="last",...){
    if (length(class(object))!=1 || class(object)!="survfit" || object$typ !="right")
        stop("Predictions only available \nfor class 'survfit', possibly stratified Kaplan-Meier fits.\n For class 'cph' Cox models see survest.cph.")
    if (missing(newdata))
        npat <- 1
    else
        if (is.data.frame(newdata))
            npat <- nrow(newdata)
        else stop("If argument `newdata' is supplied it must be a dataframe." )
    ntimes <- length(times)
    sfit <- summary(object,times=times)
    if (is.na(fill))
        Fill <- function(x,len){x[1:len]}
    else if (fill=="last")
        Fill <- function(x,len){
            y <- x[1:len]
            y[is.na(y)] <- x[length(x)]
            y}
    else stop("Argument fill must be the string 'last' or NA.")
    if (is.null(object$strata)){
        pp <- Fill(sfit$surv,ntimes)
        p <- matrix(rep(pp,npat),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    else{
        covars <- attr(terms(eval.parent(object$call$formula)),"term.labels")
        if (!all(match(covars,names(newdata),nomatch=FALSE)))
            stop("Not all strata defining variables occur in newdata.")
        ## FIXME there are different ways to build strata levels
        ## how can we test which one was used???
        stratdat <- newdata[,covars,drop=FALSE]
        names(stratdat) <- covars
        NewStratVerb <- survival::strata(stratdat)
        NewStrat <- interaction(stratdat,sep=" ")
        levs <- levels(sfit$strata)
        #    print(levs)
        #    print(levels(NewStrat))
        #    print(levels(NewStratVerb))
        if (!all(choose <- match(NewStratVerb,levs,nomatch=F))
            &&
            !all(choose <- match(NewStrat,levs,nomatch=F)))
            stop("Not all strata levels in newdata occur in fit.")
        survlist <- split(sfit$surv,sfit$strata)
        pp <- lapply(survlist[choose],Fill,ntimes)
        p <- matrix(unlist(pp,use.names=FALSE),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    if (type=="cuminc") 1-p else p
}

##' @export
predictRisk.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,type="cuminc",bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}



##' @export
predictRisk.psm <- function(object,newdata,times,...){
    if (length(times)==1){
        p <- rms::survest(object,times=c(0,times),newdata=newdata,what="survival",conf.int=FALSE)[,2]
    }else{
        p <- rms::survest(object,times=times,newdata=newdata,what="survival",conf.int=FALSE)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

##' @export 
predictRisk.rfsrc <- function(object, newdata, times, cause, ...){
    if (missing(times)||is.null(times)){
        p <- as.numeric(stats::predict(object,newdata=newdata,importance="none",...)$predicted[,2])
        ## class(p) <- "predictRisk"
        p
    }else{
        if (object$family=="surv") {
            ptemp <- predict(object,newdata=newdata,importance="none",...)$survival
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            if (!is.numeric(cause)) stop("cause is not numeric")
            cif <- predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(0,cif)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }
    }
}

##' @export 
predictRisk.FGR <- function(object,newdata,times,cause,...){
    ## require(cmprsk)
    ## predict.crr <- cmprsk:::predict.crr
    p <- predict(object=object,newdata=newdata,times=times)
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRisk.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  ## if (type=="survival")
      ## p <- cbind(1,1-temp$cuminc)[,pos+1,drop=FALSE]
  ## else
  p <- cbind(0,temp$risk)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}

##' @export 
predictRisk.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  p <- cbind(0,temp$P1)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}


##' @export 
predictRisk.CauseSpecificCox <- function (object, newdata, times, cause, se = FALSE, iid = FALSE, ...) { 
    p <- predict(object=object,newdata=newdata,times=times,cause=cause,keep.strata=FALSE, se = se, iid = iid)
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' S3-wrapper for S4 function penalized
##'
##' S3-wrapper for S4 function penalized
##' @export
##' @param formula Communicated outcome and explanatory variables. See examples.
##' @param data Data set in which formula is to be interpreted
##' @param type String specifying the type of penalization. Should match one of the following values:
##' \code{"ridge"}, \code{"lasso"}, \code{"elastic.net"}.
##' @param ... Arguments passed to penalized
##' @examples
##' \dontrun{
##' ## too slow
##' library(penalized)
##' set.seed(8)
##' d <- sampleData(200,outcome="binary")
##' newd <- sampleData(80,outcome="binary")
##' fitridge <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="ridge",
##' standardize=TRUE, model="logistic",trace=FALSE)
##' fitlasso <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="lasso",
##' standardize=TRUE, model="logistic",trace=FALSE)
##' # fitnet <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="elastic.net",
##' # standardize=TRUE, model="logistic",trace=FALSE)
##' predictRisk(fitridge,newdata=newd)
##' predictRisk(fitlasso,newdata=newd)
##' # predictRisk(fitnet,newdata=newd)
##' Score(list(fitridge),data=newd,formula=Y~1)
##' Score(list(fitridge),data=newd,formula=Y~1,splitMethod="bootcv",B=2)
##' }
##'\dontrun{ data(nki70) ## S4 fit pen <- penalized(Surv(time, event),
##' penalized = nki70[,8:77], unpenalized = ~ER+Age+Diam+N+Grade, data
##' = nki70, lambda1 = 1) penS3 <-
##' penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
##' data=nki70, lambda1=1)
##' ## or
##' penS3 <- penalizedS3(Surv(time,event)~ER+pen(TSPYL5,Contig63649_RC)+pen(10:77)+N+Grade,
##'                      data=nki70, lambda1=1)
##' ## also this works
##' penS3 <- penalizedS3(Surv(time,event)~ER+Age+pen(8:33)+Diam+pen(34:77)+N+Grade,
##'                     data=nki70, lambda1=1)
##' }
##' @export
penalizedS3 <- function(formula,data,type="ridge",...){
    # {{{ distangle the formula
    ff <- as.character(formula)
    response <- formula(paste(ff[[2]],"~1",sep=""))
    terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
    terms <- sapply(terms,function(tt){## remove whitespace
        gsub(" ","",tt)
    })
    strippedTerms <- strsplit(terms,"[()]")
    # }}}
    # {{{ extract the penalized and unpenalized parts
    penalTerms <- sapply(strippedTerms,function(x){length(x)==2 && x[[1]]=="pen"})
    unpenalVarnames <- unlist(strippedTerms[penalTerms==FALSE])
    if (length(unpenalVarnames)>0){
        unpenalized <- formula(paste("~",paste(unpenalVarnames,collapse="+")))
        response <- update.formula(response,unpenalized)
    }
    penalizedVarnames <- unlist(sapply(strippedTerms[penalTerms==TRUE],
                                       function(x){strsplit(x[[2]],",")}),use.names=FALSE)
    penalizedVarPositions <- unlist(lapply(penalizedVarnames,function(x){
        if (length(splitter <- strsplit(x,":")[[1]])>1)
            seq(as.numeric(splitter)[1],as.numeric(splitter)[2],1)
        else
            match(x,names(data),nomatch=0)
    }),use.names=FALSE)
    penalizedVarPositions <- unique(penalizedVarPositions)
    ## print(penalizedVarPositions)
    if (any(tested <- (penalizedVarPositions>NCOL(data))|penalizedVarPositions<0))
        stop("Cannot find variable(s): ",names(data[tested]))
    penalized <- data[,penalizedVarPositions]
    # }}}
    # {{{ global test
    # }}}
    type <- match.arg(tolower(type),choices=c("ridge","lasso","elastic.net"),several.ok=FALSE)
    # {{{ find optimal L1
    if (type!="ridge"){
        las1=penalized::profL1(response=response,
                    penalized=penalized,
                    data=data,
                    fold=1:nrow(data),
                    ...)
        L1=penalized::optL1(response=response,
                 penalized=penalized,
                 data=data,
                 fold=las1$fold,
                 ...)}
    # }}}
    # {{{ find optimal L2
    if (type!="lasso"){
        L2=penalized::optL2(response=response,
                 penalized=penalized,
                 data=data,
                 fold=1:nrow(data),
                 ...)}
    # }}}
    # {{{ call S4 method
    ## unpenalized terms are communicated via
    ## the left hand side of response
    fitS4 <- switch(type,"ridge"={
        penalized(response=response,
                  penalized=penalized,
                  data=data,
                  lambda1=0,
                  lambda2=L2$lambda,
                  ...)}, 
        "lasso"={
            penalized(response=response,
                      penalized=penalized,
                      data=data,
                      lambda1=L1$lambda,
                      lambda2=0,
                      ...)},
        "elastic.net"={
            penalized(response=response,
                      penalized=penalized,
                      data=data,
                      lambda1=L1$lambda,
                      lambda2=L2$lambda,
                      ...)})

    # }}}
    # {{{ deliver S3 object
    fit <- list(fitS4=fitS4,call=match.call())
    class(fit) <- "penfitS3"
    fit
    # }}}
}

##' @export 
predictRisk.penfitS3 <- function(object,
                                       newdata,
                                       ...){
    penfit <- object$fitS4
    pCovaNames <- names(penfit@penalized)
    if ("data.table" %in% class(newdata))
        newPen <- newdata[,pCovaNames,with=FALSE]
    else
        newPen <- newdata[,pCovaNames]
    p <- penalized::predict(penfit,penalized=newPen,data=newdata)
    p
}

#----------------------------------------------------------------------
### predictRisk.R ends here
