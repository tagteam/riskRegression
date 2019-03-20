### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02) 
## Version: 
## last-updated: Mar  6 2019 (18:53) 
##           By: Thomas Alexander Gerds
##     Update #: 133
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
#' predictRisk.matrix predictRisk.pecCtree
#' predictRisk.pecCforest predictRisk.prodlim predictRisk.psm
#' predictRisk.selectCox predictRisk.survfit predictRisk.randomForest
#' predictRisk.lrm predictRisk.glm
#' predictRisk.rpart
#' @usage
#' \method{predictRisk}{glm}(object,newdata,...)
#' \method{predictRisk}{cox.aalen}(object,newdata,times,...)
#' \method{predictRisk}{cph}(object,newdata,times,...)
#' \method{predictRisk}{coxph}(object,newdata,times,...)
#' \method{predictRisk}{matrix}(object,newdata,times,cause,...)
#' \method{predictRisk}{selectCox}(object,newdata,times,...)
#' \method{predictRisk}{psm}(object,newdata,times,...)
#' \method{predictRisk}{survfit}(object,newdata,times,...)
#' \method{predictRisk}{riskRegression}(object,newdata,times,cause,...)
#' \method{predictRisk}{prodlim}(object,newdata,times,cause,...)
#' \method{predictRisk}{rfsrc}(object,newdata,times,cause,...)
#' \method{predictRisk}{FGR}(object,newdata,times,cause,...)
#' \method{predictRisk}{CauseSpecificCox}(object,newdata,times,cause,...)
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return For binary outcome a vector with predicted risks. For survival outcome with and without
#' competing risks
#' a matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry is a probability and in
#' rows the values should be increasing.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @details
#' In uncensored binary outcome data there is no need to choose a time point.
#'
#' When operating on models for survival analysis (without competing risks) the function still
#' predicts the risk, as 1 - S(t|X) where S(t|X) is survival chance of a subject characterized
#' by X.
#'
#' When there are competing risks (and the data are right censored) one needs
#' to specify both the time horizon for prediction (can be a vector) and the
#' cause of the event. The function then extracts the absolute risks F_c(t|X)
#' aka the cumulative incidence of an event of type/cause c until time t for a
#' subject characterized by X. Depending on the model it may or not be possible
#' to predict the risk of all causes in a competing risks setting. For example. a
#' cause-specific Cox (CSC) object allows to predict both cases whereas a Fine-Gray regression
#' model (FGR) is specific to one of the causes. 
#' 
#' @keywords survival
#' @examples
#' ## binary outcome
#' library(rms)
#' set.seed(7)
#' d <- sampleData(80,outcome="binary")
#' nd <- sampleData(80,outcome="binary")
#' fit <- lrm(Y~X1+X8,data=d)
#' predictRisk(fit,newdata=nd)
#'\dontrun{
#' library(SuperLearner)
#' set.seed(1)
#' sl = SuperLearner(Y = d$Y, X = d[,-1], family = binomial(),
#'       SL.library = c("SL.mean", "SL.glmnet", "SL.randomForest"))
#'}
#' 
#' ## survival outcome
#' # generate survival data
##' library(prodlim)
##' set.seed(100)
##' d <- sampleData(100,outcome="survival")
##' d[,X1:=as.numeric(as.character(X1))]
##' d[,X2:=as.numeric(as.character(X2))]
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
##' ## suppose we want to predict the survival probabilities for all subjects
##' ## in the validation data at the following time points:
##' ## 0, 12, 24, 36, 48, 60
##' psurv <- predictRisk(fitCox,newdata=valdat,times=seq(0,60,12))
##' ## This is a matrix with event probabilities (1-survival)
##' ## one column for each of the 5 time points
##' ## one row for each validation set individual
##' 
##' # Do the same for a randomSurvivalForest model
##' # library(randomForestSRC)
##' # rsfmodel <- rfsrc(Surv(time,event)~X1+X2,data=learndat)
##' # prsfsurv=predictRisk(rsfmodel,newdata=valdat,times=seq(0,60,12))
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
#' train <- prodlim::SimCompRisk(100)
#' test <- prodlim::SimCompRisk(10)
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
    stop(paste0("No method available for evaluating predicted probabilities from objects in class: ",class(object),". But, you can write it yourself or ask the package manager."),call.=FALSE)
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
predictRisk.factor <- function(object,newdata,times,cause,...){
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
predictRisk.coxph <- function(object,newdata,times,...){
    p <- predictCox(object=object,
                    newdata=newdata,
                    times=times,
                    se = FALSE,
                    iid = FALSE,
                    keep.times=FALSE,
                    type="survival")$survival

    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(1-p)
}

##' @export
predictRisk.coxphTD <- function(object,newdata,times,landmark,...){
    stopifnot(attr(object$y,"type")=="counting")
    bh <- survival::basehaz(object,centered=TRUE)
    Lambda0 <- bh[,1]
    etimes <- bh[,2]
    lp <- predict(object,newdata=newdata,type="lp")
    p <- do.call("cbind",lapply(times,function(ttt){
        index <- prodlim::sindex(eval.times=c(landmark,landmark+ttt),jump.times=etimes)
        Lambda0.diff <- c(0,Lambda0)[1+index[2]] - c(0,Lambda0)[1+index[1]]
        1-exp(-Lambda0.diff * exp(lp))
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(p)
}

##' @export
predictRisk.CSCTD <- function(object,newdata,times,cause,landmark,...){
    stopifnot(attr(object$models[[1]]$y,"type")=="counting")
    if (missing(cause)) cause <- object$theCause
    else{
        ## cause <- prodlim::checkCauses(cause,object$response)
        cause <- unique(cause)
        if (!is.character(cause)) cause <- as.character(cause)
        fitted.causes <- prodlim::getStates(object$response)
        if (!(all(cause %in% fitted.causes))){
            stop(paste0("Cannot find requested cause(s) in object\n\n",
                        "Requested cause(s): ",
                        paste0(cause,collapse=", "),
                        "\n Available causes: ",
                        paste(fitted.causes,collapse=", "),"\n"))
        }
    }
    causes <- object$causes
    index.cause <- which(causes == cause)
    bh <- lapply(1:length(object$models),function(m){
        survival::basehaz(object$models[[m]],centered=TRUE)
    })
    lp <- lapply(1:length(object$models),function(m){
        predict(object$models[[m]],
                newdata=newdata,
                type="lp")
    })
    p <- do.call("cbind",lapply(times,function(ttt){
        hazard <- lapply(1:length(object$models),function(m){
            index <- sindex(eval.times=c(landmark,landmark+ttt),jump.times=bh[[m]][,2])
            bh.m <- bh[[m]][,1]
            exp(lp[[m]])*(c(0,bh.m)[1+index[2]] - c(0,bh.m)[1+index[1]])
        })
        surv <- exp(- Reduce("+",hazard))
        surv*hazard[[index.cause]]
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(p)
}


## predictRisk.coxph <- function(object,newdata,times,...){
## baselineHazard.coxph(object,times)
## require(survival)
## new feature of the survival package requires that the
## original data are included
## survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
## survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
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
predictRisk.cph <- function(object,newdata,times,...){
    ## if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
    ## p <- rms::survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
    ## if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
    p <- predictCox(object=object,
                    newdata=newdata,
                    times=times,
                    se = FALSE,
                    iid = FALSE,
                    keep.times=FALSE,
                    type="survival")$survival
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    return(1-p)
}

##' @export
predictRisk.selectCox <- function(object,newdata,times,...){
    predictRisk(object[[1]],newdata=newdata,times=times,...)
}


##' @export 
predictRisk.prodlim <- function(object,newdata,times,cause,...){
    ## require(prodlim)
    if (object$model[[1]]=="competing.risks" && missing(cause))
        stop(paste0("Cause is missing. Should be one of the following values: ",paste(attr(object$model.response,"states"),collapse=", ")))
    p <- predict(object=object,
                 cause=cause,
                 type="cuminc",
                 newdata=newdata,
                 times=times,
                 mode="matrix",
                 level.chaos=1)
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
    if (length(class(object))!=1 || class(object)[[1]]!="survfit" || object$type[[1]] !="right")
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
predictRisk.ranger <- function(object, newdata, times, cause, ...){
    if (missing(times)||is.null(times)){
        p <- stats::predict(object,data=newdata,importance="none",...)$predictions
        p
    }else{
        if (object$treetype=="Survival") {
            ptemp <- 1-stats::predict(object,data=newdata,...)$survival
            pos <- prodlim::sindex(jump.times=ranger::timepoints(object),eval.times=times)
            p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            message("Sorry. Hope Marvin does this soon.")
        }
    }
}

##' @export 
predictRisk.rfsrc <- function(object, newdata, times, cause, ...){
    if (missing(times)||is.null(times)){
        p <- as.numeric(stats::predict(object,newdata=newdata,importance="none",...)$predicted[,2])
        ## class(p) <- "predictRisk"
        p
    }else{
        if (object$family=="surv") {
            ptemp <- 1-stats::predict(object,newdata=newdata,importance="none",...)$survival
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            if (is.character(cause)) cause <- as.numeric(cause)
            if (!is.numeric(cause)) stop("cause is not numeric")
            cif <- stats::predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
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
predictRisk.CauseSpecificCox <- function (object, newdata, times, cause, ...) { 
    p <- predict(object=object,
                 newdata=newdata,
                 times=times,
                 cause=cause,
                 keep.strata=FALSE,
                 se = FALSE,
                 iid = FALSE)$absRisk
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export 
predictRisk.BinaryTree<-function (object, newdata, times, ...) 
{
    learndat <- mod@data@get("input")
    learnresponse <- mod@data@get("response")
    resp <- as.data.frame(learnresponse)
    browser()
    nclass <- length(unique(object@where))
    learndat$ctreeFactor <- factor(predict(object, newdata = learndat,type="node")) #was vector
    newdata$ctreeFactor <- factor(predict(object, newdata = newdata,type="node")) #was vector
    BinaryTree.form <- reformulate("ctreeFactor",object@data@formula$response[[2]])
    fit.BinaryTree <- prodlim(BinaryTree.form, data = learndat)
    p <- predictRisk(fit.BinaryTree, newdata = newdata, times = times) #input prodlim object
    1 - p
}


##' @export 
predict.smcure <- function(object, newdata, times,  ...) {
    ##Wrapper for predictsmcure so the X's and Z's do not need to be prespecified.
    require(dummies)
    # pred <- predictsmcure(object, newX = dummy.data.frame(newdata)[,object$betanm], newZ = dummy.data.frame(newdata)[,object$bnm[-1]], model = "ph")$prediction  ##object$call$model -> werkt niet!!
    pred <- predictsmcure(object, newX = dummy.data.frame(newdata)[,object$betanm[-1]], newZ = dummy.data.frame(newdata)[,object$bnm[-1]], model = "ph")$prediction  ##object$call$model -> werkt niet!!
    pdsort <- pred[order(pred[, "Time"]), ]
    pdsort <- rbind(c(array(1,ncol(pdsort)-1),0),pdsort) #rij met tijd=0 en surv = 1 toevoegen
    orgTimes <- pdsort[,"Time"]
    tim <- array(0,length(times))       #index of time-points
    for (i in 1:length(times)){         #reproduce step function of S0 estimation sample
        tim[i] <- max(which(orgTimes-times[i]<=0))}
    p <- t(pdsort[tim,-ncol(pdsort)])   #leave out last column (time)
    1 - p
}

##' @export 
predictRisk.smcure <- function (object, newdata, times, ...) {
  p <- predict.smcure(object, newdata = newdata, times = times )
}

##' @export 
predictRisk.penfit <- function(object, newdata, times, ...){
    ##only works with a optL1 or optL2 object. 
    require(penalized)
    pred1 <- predict(object, penalized=object@formula$penalized, data = newdata)
    p <- matrix(0,nrow(newdata),length(times))
    for (i in 1:length(times)){
        p[,i] <- survival(pred1, time = times[i])
    }
    1 - p
}

##' @export 
predictRisk.flexsurvreg <- function(object, newdata, times, ...) {
    require(dummies)
    require(flexsurv)
    p <- matrix(0, NROW(newdata), length(times))
    sm <- summary.flexsurvreg(object, X = as.matrix(dummy.data.frame(newdata)[,dimnames(object$data$X)[[2]]]), t = times, start = 0, B = 0) #no confidence interval simulations
    for (i in 1:NROW(newdata)){
        p[i,] <- sm[[i]][,2]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop("Prediction failed")
    1 - p
}

scaleorg <- function(Xt, X){
    ##scale new X matrix like training matrix.
    ##Xt is X_train, X is new X
    xtmean <- apply(Xt, MARGIN = 2, mean ) #mean of training matrix
    Xscaled <- sweep(X, 2L, xtmean, check.margin = FALSE)
    f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/max(1, length(v) - 1L))
    }
    xtsd <- apply(Xt, 2L, f)             #sd of training matrix
    Xscaled <- sweep(X, 2L, xtsd, "/", check.margin = FALSE)
    return(Xscaled)
}

predict.coxpls2 <- function(object, newdata, traindata, ...){
  require(dummies)
  vars <- dimnames(object$pls2$pls_mod$coefficients)[[1]]
  loading.weights <- object$pls2$pls_mod$loading.weights
  if (identical(newdata, traindata)){
    scores <- scale(data.matrix(dummy.data.frame(traindata))[, vars]) %*% loading.weights
  }
  else {
    Xscaled <- scaleorg(data.matrix(dummy.data.frame(traindata))[, vars], data.matrix(dummy.data.frame(newdata))[, vars])
    scores <- Xscaled %*% loading.weights
  }
  colnames(scores) <- paste("Comp",seq(1,object$pls2$pls_mod$ncomp,1),sep=".")
  sf <- survfit(object$pls2$cox_pls, newdata = as.data.frame(scores))
  ## sf <- survfit(object$cox_pls, newdata = scores)
  sf <- cbind(sf$time,sf$surv)
  return(sf)
}

plscox <- function (formula, data, ncomp = 3, cv = FALSE, ...) 
{
    ## formula interface to coxpls2
    require(plsRcox)
    call <- match.call(expand.dots = FALSE)
    formula.names <- try(all.names(formula), silent = TRUE)
    actual.terms <- terms(formula, data = data)
    formula <- eval(call$formula)
    response <- model.response(model.frame(formula, data))
    Time <- as.numeric(response[, "time"])
    Event <- as.numeric(response[, "status"])
    X <- model.matrix(actual.terms, data = data)[, -c(1), drop = FALSE]
    nc <<- ncomp                    #to prevent namespace problems
    pl <- coxpls2(X, time = Time, time2 = Event, allres = TRUE, ncomp = nc)
    out <- list(pls2 = pl, ncomp = ncomp, 
                call = call, formula = formula, response = response, vars = attr(actual.terms,"term.labels"))
    class(out) <- "coxpls2"
    return(out)
}

##' @export 
predictRisk.coxpls2 <- function(object, newdata, times, traindata = NULL){
    require(plsRcox)
    ## calculate componentscores on basis of newdata
    ## apply survfit with newdata = Comp.1 t/m Comp.k
    pdsort <- predict.coxpls2(object, newdata = newdata, traindata)
    pdsort <- rbind(c(0,array(1,ncol(pdsort)-1)),pdsort) #add row with time =0 and surv = 1 
    orgTimes <- pdsort[,1]                               #time in 1st column
    tim <- array(0,length(times))       #index of time-points
    for (i in 1:length(times)){         #reproduce step function of S0 of estimation sample 
        tim[i] <- max(which(orgTimes-times[i]<=0))}
    p <- t(pdsort[tim,-1])
    1 - p
}

##' @export 
predictRisk.gbm <- function(object, newdata, times, traindata, n.trees = NULL) {
    ## n.trees passed via ...?
    if (is.null(n.trees)) n.trees <- object$n.trees
    p <- matrix(0, NROW(newdata), length(times))
    xb.train <- predict(object ,newdata = traindata, n.trees = n.trees)
    H2 <- basehaz.gbm(t = traindata[,object$response.name[2]],delta = traindata[,object$response.name[3]], f.x = xb.train, t.eval = times)
    xb.test <- predict(object, newdata = newdata , n.trees = n.trees ) 
    for (i in 1:length(times)) p[,i] <- exp(-H2[i] * exp(xb.test))
    p[,times==0] <- 1                   #to prevent problems with prediction error curves 
    1 - p
}

##' @export 
predictRisk.phnnet <- function (object, newdata, times, traindata = NULL){
    require(survnnet)
    learndat <- object$traindata
    seeds <- sample(1:1000, size = 3)
    object$call$data <- learndat
    re.fitter <- lapply(seeds, function(s) {
        set.seed(s)
        refit <- eval(object$call)
        list(learn = predict(refit, learndat), val = predict(refit, 
                                                             newdata))
    })
    learndat$nnetFactor <- rowMeans(do.call("cbind", lapply(re.fitter, 
                                                            function(x) x[["learn"]])))
    newdata$nnetFactor <- rowMeans(do.call("cbind", lapply(re.fitter, 
                                                           function(x) x[["val"]])))
    nnet.form <- reformulate("nnetFactor", object$call$formula[[2]])
    fit.nnet <- cph(nnet.form, data = learndat, se.fit = FALSE, 
                    surv = TRUE, x = TRUE, y = TRUE)
    p <- predictRisk.cph(fit.nnet, newdata = newdata, times = times)
    1 - p
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
##' library(prodlim)
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
##' Score(list(fitridge),data=newd,formula=Y~1,split.method="bootcv",B=2)
##' }
##'\dontrun{ data(nki70) ## S4 fit
##' pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
##'                  unpenalized = ~ER+Age+Diam+N+Grade, data = nki70,
##' lambda1 = 1)
##' penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
##'                      data=nki70, lambda1=1)
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


##' @export 
predictRisk.SuperPredictor  <- function(object,newdata,...){
    p <- SuperLearner::predict.SuperLearner(object=object,newdata=newdata)$pred
    ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        ## stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}



#----------------------------------------------------------------------
### predictRisk.R ends here
