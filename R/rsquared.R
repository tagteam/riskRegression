### rsquared.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug  9 2017 (10:36) 
## Version: 
## Last-Updated: Aug 10 2017 (18:36) 
##           By: Thomas Alexander Gerds
##     Update #: 62
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
##' @param object Model for which we want R^2 
##' @param ... passed to \code{riskRegression::Score}
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
rsquared.default <- function(object,times,...){
  stop("No method available for calculating R^2 for objects in class: ",class(object),call.=FALSE)
}

##' @export
rsquared.CauseSpecificCox <- function(object,times,...){
    Rsquared=model=Variable=Rsquared.gain=Brier=NULL
    ## this is a modification of drop1
    models <- object$models
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
    r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,Rsquared,Rsquared.gain)])
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
    r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,Rsquared,Rsquared.gain)])
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
    r2[,Rsquared:=100*Rsquared]
    ## r2 <- r2[model!="Null model"]
    data.table(setnames(r2,"model","Variable"))
    r2[,Rsquared.gain:=Rsquared[Variable=="Full model"]-Rsquared]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,Brier,Rsquared,Rsquared.gain)])
    out
}


######################################################################
### rsquared.R ends here
