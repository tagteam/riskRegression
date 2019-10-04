### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02) 
## Version: 
## last-updated: okt  4 2019 (15:58) 
##           By: Brice Ozenne
##     Update #: 264
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
#' predictRisk.rpart predictRisk.gbm
#' predictRisk.flexsurvreg
#' @usage
#' \method{predictRisk}{glm}(object,newdata,...)
#' \method{predictRisk}{cox.aalen}(object,newdata,times,...)
#' \method{predictRisk}{cph}(object,newdata,times,product.limit,...)
#' \method{predictRisk}{coxph}(object,newdata,times,product.limit,...)
#' \method{predictRisk}{matrix}(object,newdata,times,cause,...)
#' \method{predictRisk}{selectCox}(object,newdata,times,...)
#' \method{predictRisk}{psm}(object,newdata,times,...)
#' \method{predictRisk}{survfit}(object,newdata,times,...)
#' \method{predictRisk}{riskRegression}(object,newdata,times,cause,...)
#' \method{predictRisk}{prodlim}(object,newdata,times,cause,...)
#' \method{predictRisk}{rfsrc}(object,newdata,times,cause,...)
#' \method{predictRisk}{FGR}(object,newdata,times,cause,...)
#' \method{predictRisk}{CauseSpecificCox}(object,newdata,times,cause,product.limit,...)
#' \method{predictRisk}{gbm}(object,newdata,times,n.trees,...)
#' \method{predictRisk}{flexsurvreg}(object,newdata,times,...)
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param iid Should the iid decomposition be output using an attribute?
#' @param average.iid Should the average iid decomposition be output using an attribute?
#' @param product.limit If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
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
predictRisk.glm <- function(object, newdata, iid = FALSE, average.iid = FALSE,...){

    if (object$family$family=="binomial"){
        if(iid & average.iid){
            stop("Cannot specify both \'iid\'=TRUE and \'average.iid\'=TRUE \n")
        }else if(iid){
            out <- predictGLM(object,
                              newdata = newdata,
                              average.iid = average.iid)
        }else if(average.iid){
            out <- predictGLM(object,
                              newdata = newdata,
                              average.iid = average.iid)
        }else{
            out <- as.numeric(stats::predict(object,newdata=newdata,type="response"))
        }
        ## hidden argument: enable to ask for the prediction of Y==1 or Y==0
        level <- list(...)$level
        if(!is.null(level)){
            matching.Ylevel <- table(object$data[[all.vars(formula(object))[1]]],
                                     object$y)
            all.levels <- rownames(matching.Ylevel)
            level <- match.arg(level, all.levels)

            index.level <- which(matching.Ylevel[level,]>0)
            if(length(index.level) > 1){
                stop("Unknown value for the outcome variable \n")
            }else if(index.level == 1){
                out <- 1 - out
                if(iid){
                    attr(out,"iid") <- - attr(out,"iid")
                }else if(average.iid){
                    if(is.list(attr(out,"iid"))){
                        attr(out,"average.iid") <- lapply(attr(out,"iid"), function(iIID){-iIID})
                    }else{
                        attr(out,"average.iid") <- - attr(out,"iid")
                    }
                }                
            }else if(average.iid){
                attr(out,"average.iid") <- attr(out,"iid")
            }
        }else if(average.iid){
            attr(out,"average.iid") <- attr(out,"iid")
        }
        return(out)
    } else {
        stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
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
predictRisk.coxph <- function(object, newdata, times, product.limit = FALSE, iid = FALSE, average.iid = FALSE, ...){
    type <- list(...)$type ## hidden argument for ate
    
    if(product.limit){
        outPred <- predictCoxPL(object=object,
                                newdata=newdata,
                                times=times,
                                iid = iid,
                                iid = average.iid,
                                keep.times=FALSE,
                                type="survival")
    }else{
        outPred <- predictCox(object=object,
                              newdata=newdata,
                              times=times,
                              iid = iid,
                              average.iid = average.iid,
                              keep.times=FALSE,
                              type="survival")
    }
    if(identical(type,"survival")){
        out <- outPred$survival
    }else{
        out <- 1-outPred$survival
    }
    if(iid){
        if(identical(type,"survival")){
            attr(out,"iid") <- outPred$survival.iid
        }else{
            attr(out,"iid") <- -outPred$survival.iid
        }
    }
    if(average.iid){
        if(identical(type,"survival")){
            attr(out,"average.iid") <- outPred$survival.average.iid
        }else{
            if(is.list(outPred$survival.average.iid)){
                attr(out,"average.iid") <- lapply(outPred$survival.average.iid, function(iIID){-iIID})
            }else{
                attr(out,"average.iid") <- -outPred$survival.average.iid
            }
        }
    }

    return(out)
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
predictRisk.cph <- predictRisk.coxph

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
    xvars <- object$forest$independent.variable.names
    newdata <- subset(newdata,select=xvars)
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
predictRisk.CauseSpecificCox <- function (object, newdata, times, cause, product.limit = TRUE, iid = FALSE, average.iid = FALSE, ...) { 
    diag <- list(...)$diag
    if(is.null(diag)){diag <- FALSE}
    
    outPred <- predict(object=object,
                       newdata=newdata,
                       times=times,
                       cause=cause,
                       keep.strata=FALSE,
                       se = FALSE,
                       iid = iid,
                       diag = diag,
                       average.iid = average.iid,
                       product.limit = product.limit)
    
    out <- outPred$absRisk
    if(iid){
        attr(out,"iid") <- outPred$absRisk.iid
    }
    if(average.iid){
        attr(out,"average.iid") <- outPred$absRisk.average.iid
    }

    return(out)
}




##' S3-wrapper for S4 function penalized
##'
##' S3-wrapper for S4 function penalized
##' @export
##' @param formula Communicated outcome and explanatory variables. See examples.
##' @param data Data set in which formula is to be interpreted
##' @param type String specifying the type of penalization. Should match one of the following values:
##' \code{"ridge"}, \code{"lasso"}, \code{"elastic.net"}.
##' @param lambda1 Lasso penalty
##' @param lambda2 ridge penalty
##' @param fold passed to \code{penalized::profL1}
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
penalizedS3 <- function(formula,
                        data,
                        type="elastic.net",
                        lambda1,
                        lambda2,
                        fold,
                        ...){
                                        # {{{ distangle the formula
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=c("pen","unpen"),
                                       stripSpecials=c("pen","unpen"),
                                       stripUnspecials="pen",
                                       specialsDesign=TRUE)
    response <- EHF$event.history
    pen <- EHF$pen
    unpen <- EHF$unpen
    if (is.null(unpen))
        args <- list(response=response,penalized=pen,data=data,...)
    else
        args <- list(response=response,penalized=pen,unpenalized=unpen,data=data,...)
                                        # {{{ find optimal L1
    if ((tolower(type) %in% c("lasso","elastic.net")) && missing(lambda1)){
        if (missing(fold)) fold <- 1:NROW(data)
        las1 <- do.call(penalized::profL1,c(args,list(fold=fold)))
        lambda1 <- do.call(penalized::optL1,c(args,list(fold=las1$fold)))$lambda
    }
                                        # }}}
                                        # {{{ find optimal L2
    if ((tolower(type) %in% c("ridge","elastic.net")) && missing(lambda2)){
        if (missing(fold)) fold <- 1:NROW(data)
        lambda2 <- do.call(penalized::optL2,c(args,list(fold=fold)))$lambda
    }
                                        # }}}
                                        # {{{ call S4 method
    ## unpenalized terms are communicated via
    ## the left hand side of response
    fitS4 <- switch(type,"ridge"={
        do.call(penalized::penalized,c(args,list(lambda1=0, lambda2=lambda2)))
    }, 
    "lasso"={
        do.call(penalized::penalized,c(args,list(lambda1=lambda1, lambda2=0)))
    },
    "elastic.net"={
        do.call(penalized::penalized,c(args,list(lambda1=lambda1, lambda2=lambda2)))
    })
    if (is.infinite(fitS4@loglik)){
        print(fitS4)
        stop()
    }
                                        # }}}
                                        # {{{ deliver S3 object
    fit <- list(fitS4=fitS4,call=match.call())
    fit$terms <- terms(formula)
    class(fit) <- "penfitS3"
    fit
                                        # }}}
}

##' @export 
predictRisk.penfitS3 <- function(object,
                                 newdata,
                                 times,
                                 ...){
    penfit <- object$fitS4
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (NROW(newdata) == 0) stop("No (non-missing) observations")
    rhs <- as.formula(delete.response(object$terms))
    newdata$dummy.time=1
    newdata$dummy.event=1
    dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    EHF <- prodlim::EventHistory.frame(dummy.formula,
                                       data=newdata,
                                       specials=c("pen","unpen"),
                                       stripSpecials=c("pen","unpen"),
                                       stripUnspecials="pen",
                                       specialsDesign=TRUE)
    pen <- EHF$pen
    unpen <- EHF$unpen
    args <- list(penfit,penalized=pen)
    if (length(unpen)>0) args <- c(args,list(unpenalized=unpen))
    p <- do.call(penalized::predict,args)
    if (penfit@model=="cox"){
        if (missing(times)) stop("Need time points for predicting risks from a penalized Cox regression model.")
        ttt <- p@time
        p <- cbind(0,1-p@curves)[,prodlim::sindex(jump.times=ttt,eval.times=times)+1,drop=FALSE]
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        }
    }
    p
}

##' @title SmcFcs 
##' @description TODO
##' 
##' @param formula TODO
##' @param data TODO
##' @param m TODO
##' @param method TODO
##' @param fitter TODO
##' @param fit.formula TODO
##' @param ... TODO
##' 
##' @export
SmcFcs  <- function(formula,data,m=5,method,fitter="glm",fit.formula,...){
    requireNamespace("smcfcs")
    this <- as.character(formula)
    Yname <- this[2]
    Xnames <- this[3]
    sform <- paste0(Yname,"~",Xnames)
    if (length(unique(data[[Yname]]))!=2)
        stop("Outcome must be binary")
    Xvars <- all.vars(formula)
    if (missing(fit.formula)) fit.formula <- formula
    Avars <- all.vars(fit.formula)
    Xvars <- union(Xvars,Avars)
    Xdata <- subset(data,select=Xvars)
    if (missing(method)){
        method <- sapply(Xdata,function(x){
            if (any(is.na(x))){
                if(length(unique(x))==2){
                    "logreg"
                }else{
                    if(is.factor(x)){
                        "mlogit"
                    } else{
                        "norm"
                    }
                }
            } else{
                ""
            }
        })
    }
    idata.list <- smcfcs::smcfcs(smformula=sform,
                                 originaldata=Xdata,
                                 m=m,
                                 smtype="logistic",...,
                                 method=method)$impDatasets
    res <- lapply(idata.list,function(d){
        do.call(fitter,list(fit.formula,data=d,family="binomial"))
    })
    class(res) <- "SmcFcs"
    res
}

##' @export
predictRisk.SmcFcs <- function(object,newdata,...){
    p <- Reduce("+",lapply(object,function(x){
        predictRisk(x,newdata=newdata)
    }))/length(object)
    p
}

##' @export 
predictRisk.SuperPredictor  <- function(object,newdata,...){
    p <- SuperLearner::predict.SuperLearner(object=object,newdata=newdata)$pred
    ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        ## stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export 
predictRisk.gbm <- function(object, newdata, times, ...) {
    n.trees <- object$n.trees
    traindata <-  gbm::reconstructGBMdata(object)
    p <- matrix(0, NROW(newdata), length(times))
    xb.train <- predict(object ,newdata = traindata, n.trees = n.trees)
    H2 <- gbm::basehaz.gbm(t = traindata[, as.character(object$call$formula[[2]][[2]])], 
                      delta = traindata[, as.character(object$call$formula[[2]][[3]])], 
                      f.x = xb.train, t.eval = times)
    xb.test <- predict(object, newdata = newdata , n.trees = n.trees ) 
    for (i in 1:length(times)) p[,i] <- exp(-H2[i] * exp(xb.test))
    p[,times==0] <- 1
    1 - p
}
##' @export 
predictRisk.flexsurvreg <- function(object, newdata, times, ...) {
    newdata <- data.frame(newdata)
    p <- matrix(0, NROW(newdata), length(times))
    term <- attr(terms(as.formula(object$call$formula)), "term.labels")
    sm <- summary(object, newdata = newdata[, term], t = times, start = 0, B = 0) #no confidence interval simulations
    for (i in 1:NROW(newdata)){
        p[i,] <- sm[[i]][,2]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop("Prediction failed")
    1 - p
}

#----------------------------------------------------------------------
### predictRisk.R ends here

