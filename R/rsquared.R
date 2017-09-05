### rsquared.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug  9 2017 (10:36) 
## Version: 
## Last-Updated: Sep  5 2017 (08:56) 
##           By: Thomas Alexander Gerds
##     Update #: 95
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' General R^2 for binary outcome and right censored time to event (survival) outcome also with competing risks
##'
##' R^2 is calculated based on the model's predicted risks. The Brier score of the model is compared to the Brier score of the null model.
##' @title Explained variation for settings with binary, survival and competing risk outcome
##' @aliases rsquared rsquared.default rsquared.glm rsquared.coxph rsquared.CauseSpecificCox
##' @param object Model for which we want R^2
##' @param cause For competing risk models the event of interest
##' @param times Vector of time points used as prediction horizon for the computation of Brier scores. 
##' @param ... passed to \code{riskRegression::Score}
##' @usage
##' rsquared(object,...)
##' \method{rsquared}{default}(object,...)
##' \method{rsquared}{glm}(object,...)
##' \method{rsquared}{coxph}(object,times,...)
##' \method{rsquared}{CauseSpecificCox}(object,times,cause,...)
##' 
##' @return Data frame with explained variation values for the full model.
##' @seealso Score
##' @examples
##'
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(48,outcome="binary")
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' rsquared(lr1)
##'
##' library(survival)
##' data(pbc)
##' cox1= coxph(Surv(time,status!=0)~age+sex+log(bili)+log(albumin)+log(protime),
##' data=na.omit(pbc),x=TRUE)
##' rsquared(cox1,times=1000)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
rsquared <- function(object,...){
    UseMethod("rsquared",object)
}

##' @export
rsquared.default <- function(object,...){
  stop("No method available for calculating R^2 for objects in class: ",class(object),call.=FALSE)
}

##' @export
rsquared.CauseSpecificCox <- function(object,times,cause,...){
    Rsquared=model=Variable=Rsquared.gain=Brier=NULL
    ## this is a modification of drop1
    models <- object$models
    Terms <- lapply(models,stats::terms)
    ## tl <- attr(Terms, "term.labels")
    scope.list <- lapply(models,stats::drop.scope)
    scope <- unique(unlist(scope.list,use.names=0))
    ns <- length(scope)
    n0 <- lapply(models,stats::nobs, use.fallback = TRUE)
    leaveOneOut <- lapply(scope,function(var){
        varfit <- object
        varfit$models <- lapply(1:length(models),function(m){
            fit.m <- stats::update(models[[m]], as.formula(paste("~ . -", var)), 
                                   evaluate = FALSE)
            env <- environment(formula(models[[m]]))
            fit.m <- eval(fit.m, envir = env)
            nnew <- stats::nobs(models[[m]], use.fallback = TRUE)
            if (all(is.finite(c(n0[[m]], nnew))) && nnew != n0[[m]])
                stop("number of rows in use has changed: remove missing values?")
            fit.m
        })
        varfit
    })
    names(leaveOneOut) <- scope
    ## fixme: could constract a formula with all variables from all causes
    if (missing(cause)) cause <- attr(object$response,"states")[[1]]
    else cause <- prodlim::checkCauses(cause,object$response)
    r2 <- Score(c("Full model"=list(object),leaveOneOut),
                formula=eval(object$call$formula)[[1]],
                data=eval(object$call$data),
                times=times,
                cause=cause,
                contrasts=FALSE,
                conf.int=FALSE,
                metrics="brier",
                censModel="km",
                summary="rsquared",
                ...)$Brier$score
    ## r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,Rsquared,Rsquared.gain)])
    class(out) <- c("rsquared",class(out))
    out
}

##' @export
rsquared.coxph <- function(object,times,...){
    Rsquared=model=Variable=Rsquared.gain=Brier=NULL
    ## this is a modification of drop1
    Terms <- stats::terms(object)
    ## tl <- attr(Terms, "term.labels")
    scope <- stats::drop.scope(object)
    ns <- length(scope)
    n0 <- stats::nobs(object, use.fallback = TRUE)
    env <- environment(formula(object))
    leaveOneOut <- lapply(scope,function(var){
        varfit <- stats::update(object, as.formula(paste("~ . -", var)), 
                         evaluate = FALSE)
        varfit <- eval(varfit, envir = env)
        nnew <- stats::nobs(varfit, use.fallback = TRUE)
        if (all(is.finite(c(n0, nnew))) && nnew != n0)
            stop("number of rows in use has changed: remove missing values?")
        varfit
    })
    names(leaveOneOut) <- scope
    wdat <- eval(object$call$data)
    ## r2 <- Score(leaveOneOut,formula=Surv(time,status!=1)~1,data=pbc,times=times,metrics="rsquared")
    r2 <- Score(c("Full model"=list(object),leaveOneOut),formula=formula(object),data=wdat,times=times,contrasts=FALSE,conf.int=FALSE,metrics="brier",censModel="km",summary="rsquared",...)$Brier$score
    ## r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,Rsquared,Rsquared.gain)])
    class(out) <- c("rsquared",class(out))
    out
}

##' @export
rsquared.glm <- function(object,...){
    Rsquared=model=Variable=Rsquared.gain=Brier=NULL
    ## this is a modification of drop1
    stopifnot(family(object)$family=="binomial")
    Terms <- stats::terms(object)
    ## tl <- attr(Terms, "term.labels")
    scope <- stats::drop.scope(object)
    ns <- length(scope)
    n0 <- stats::nobs(object, use.fallback = TRUE)
    env <- environment(formula(object))
    leaveOneOut <- lapply(scope,function(var){
        varfit <- stats::update(object, as.formula(paste("~ . -", var)), 
                                evaluate = FALSE)
        varfit <- eval(varfit, envir = env)
        nnew <- stats::nobs(varfit, use.fallback = TRUE)
        if (all(is.finite(c(n0, nnew))) && nnew != n0)
            stop("number of rows in use has changed: remove missing values?")
        varfit
    })
    names(leaveOneOut) <- scope
    wdat <- eval(object$call$data)
    r2 <- Score(c("Full model"=list(object),leaveOneOut),formula=formula(object),data=wdat,contrasts=FALSE,conf.int=FALSE,metrics="brier",summary="rsquared",...)$Brier$score
    ## r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,Brier,Rsquared,Rsquared.gain)])
    class(out) <- c("rsquared",class(out))
    out
}


######################################################################
### rsquared.R ends here
