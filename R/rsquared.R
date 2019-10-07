### rsquared.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Aug  9 2017 (10:36) 
## Version: 
## Last-Updated: okt  7 2019 (18:58) 
##           By: Brice Ozenne
##     Update #: 153
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
##' @aliases rsquared rsquared.default rsquared.glm rsquared.coxph rsquared.CauseSpecificCox IPA IPA.default IPA.glm IPA.coxph IPA.CauseSpecificCox
##' @param object Model for which we want R^2
##' @param newdata Optional validation data set in which to compute R^2
##' @param formula Formula passed to \code{Score}. If not provided, try to use the formula of the call of \code{object}, if any.
##' @param cause For competing risk models the event of interest
##' @param times Vector of time points used as prediction horizon for the computation of Brier scores.
##' @param ... passed to \code{riskRegression::Score}
##' @usage
##' rsquared(object,...)
##' IPA(object,...)
##' \method{rsquared}{default}(object,formula,newdata,times,cause,...)
##' \method{rsquared}{glm}(object,formula,newdata,...)
##' \method{rsquared}{coxph}(object,formula,newdata,times,...)
##' \method{rsquared}{CauseSpecificCox}(object,formula,newdata,times,cause,...)
##' \method{IPA}{default}(object,formula,newdata,times,cause,...)
##' \method{IPA}{glm}(object,formula,newdata,...)
##' \method{IPA}{coxph}(object,formula,newdata,times,...)
##' \method{IPA}{CauseSpecificCox}(object,formula,newdata,times,cause,...)
##' 
##' @return Data frame with explained variation values for the full model.
##' @seealso Score
##' @examples
##' library(prodlim)
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(48,outcome="binary")
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' rsquared(lr1)
##' 
##' ## validation data
##' valdat=sampleData(94,outcome="binary")
##' rsquared(lr1,newdata=valdat)
##'
##' ## predicted risks externally given
##' p1=predictRisk(lr1,newdata=valdat)
##' rsquared(p1,formula=Y~1,valdat)
##' 
##' # survival
##' library(survival)
##' data(pbc)
##' pbc=na.omit(pbc)
##' pbctest=(1:NROW(pbc)) %in% sample(1:NROW(pbc),size=.632*NROW(pbc))
##' pbclearn=pbc[pbctest,]
##' cox1= coxph(Surv(time,status!=0)~age+sex+log(bili)+log(albumin)+log(protime),
##'       data=pbclearn,x=TRUE)
##' 
##' ## same data
##' rsquared(cox1,formula=Surv(time,status!=0)~1,times=1000)
##'
##' ## validation data
##' pbcval=pbc[!pbctest,]
##' rsquared(cox1,formula=Surv(time,status!=0)~1,newdata=pbcval,times=1000)
##'
##' ## predicted risks externally given
##' p2=predictRisk(cox1,newdata=pbcval,times=1000)
##' rsquared(cox1,formula=Surv(time,status!=0)~1,newdata=pbcval,times=1000)
##'  
##' # competing risks
##' data(Melanoma)
##' Melanomatest=(1:NROW(Melanoma)) %in% sample(1:NROW(Melanoma),size=.632*NROW(Melanoma))
##' Melanomalearn=Melanoma[Melanomatest,]
##' fit1 <- CSC(list(Hist(time,status)~sex,
##'                  Hist(time,status)~invasion+epicel+age),
##'                  data=Melanoma)
##' rsquared(fit1,times=1000,cause=2)
##'
##' ## validation data
##' Melanomaval=Melanoma[!Melanomatest,]
##' rsquared(fit1,formula=Hist(time,status)~1,newdata=Melanomaval,times=1000)
##'
##' ## predicted risks externally given
##' p3= predictRisk(fit1,cause=1,newdata=Melanomaval,times=1000)
##' rsquared(p3,formula=Hist(time,status)~1,cause=1,newdata=Melanomaval,times=1000)
##'  
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
rsquared <- function(object,...){
    UseMethod("rsquared")
}

##' @export
rsquared.default <- function(object,formula,newdata,times,cause,...){
    if (missing(formula)) stop("Need formula to define the outcome variable.")
    if (missing(newdata)) stop("Need newdata to calculate R^2.")
    ## stop("No method available for calculating R^2 for objects in class: ",class(object),call.=FALSE)
    if (missing(times)) times=NULL
    if (missing(cause)) cause=NULL
    out <- Score(list(object),
                 formula=formula,
                 data=newdata,
                 times=times,
                 cause=cause,
                 contrasts=FALSE,
                 conf.int=FALSE,
                 metrics="brier",
                 cens.model="km",
                 summary="rsquared",
                 ...)$Brier$score
    class(out) <- c("IPA",class(out))
    out
}

##' @export
rsquared.CauseSpecificCox <- function(object,formula,newdata,times,cause,...){
    IPA=model=Variable=IPA.drop=Brier=NULL
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
            env <- environment(models[[m]]$formula)
            fit.m <- eval(fit.m, envir = env)
            nnew <- stats::nobs(models[[m]], use.fallback = TRUE)
            if (all(is.finite(c(n0[[m]], nnew))) && nnew[[1]] != n0[[m]])
                stop("number of rows in use has changed: remove missing values?")
            fit.m
        })
        varfit
    })
    names(leaveOneOut) <- scope
    ## fixme: could constract a formula with all variables from all causes
    if (missing(cause)) cause <- attr(object$response,"states")[[1]]
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
    if (missing(newdata)) newdata=eval(object$call$data)
    if (missing(formula)) formula=eval(object$call$formula)
    ## CSC may have two formulas
    if (is.list(formula)) formula=formula[[1]]
    r2 <- Score(c("Full model"=list(object),leaveOneOut),
                formula=formula,
                data=newdata,
                times=times,
                cause=cause,
                contrasts=FALSE,
                conf.int=FALSE,
                metrics="brier",
                cens.model="km",
                summary="rsquared",
                ...)$Brier$score
    ## r2[,IPA:=100*IPA]
    ## r2 <- r2[model!="Null model"]
    data.table::setnames(r2,"model","Variable")
    r2[,IPA.drop:=IPA[Variable=="Full model"]-IPA,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,IPA,IPA.drop)])
    class(out) <- c("IPA",class(out))
    out
}

##' @export
rsquared.coxph <- function(object,formula,newdata,times,...){
    IPA=model=Variable=IPA.drop=Brier=NULL
    ## this is a modification of drop1
    Terms <- stats::terms(object)
    ## tl <- attr(Terms, "term.labels")
    scope <- stats::drop.scope(object)
    ns <- length(scope)
    n0 <- stats::nobs(object, use.fallback = TRUE)
    env <- environment(object$formula)
    leaveOneOut <- lapply(scope,function(var){
        varfit <- stats::update(object, as.formula(paste("~ . -", var)), 
                         evaluate = FALSE)
        varfit <- eval(varfit, envir = env)
        nnew <- stats::nobs(varfit, use.fallback = TRUE)
        if (all(is.finite(c(n0, nnew))) && nnew[[1]] != n0)
            stop("number of rows in use has changed: remove missing values?")
        varfit
    })
    names(leaveOneOut) <- scope
    if (missing(newdata)) newdata=eval(object$call$data)
    if (missing(formula)) formula=object$formula
    r2 <- Score(c("Full model"=list(object),leaveOneOut),
                formula=formula,
                data=newdata,
                times=times,
                contrasts=FALSE,
                conf.int=FALSE,
                metrics="brier",
                cens.model="km",
                summary="rsquared",
                ...)$Brier$score
    ## r2[,IPA:=100*IPA]
    ## r2 <- r2[model!="Null model"]
    data.table::setnames(r2,"model","Variable")
    r2[,IPA.drop:=IPA[Variable=="Full model"]-IPA,by=times]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,times,Brier,IPA,IPA.drop)])
    class(out) <- c("IPA",class(out))
    out
}

##' @export
rsquared.glm <- function(object,formula,newdata,...){
    IPA=model=Variable=IPA.drop=Brier=NULL
    ## this is a modification of drop1
    stopifnot(family(object)$family=="binomial")
    Terms <- stats::terms(object)
    ## tl <- attr(Terms, "term.labels")
    scope <- stats::drop.scope(object)
    ns <- length(scope)
    n0 <- stats::nobs(object, use.fallback = TRUE)
    env <- environment(object$formula)
    leaveOneOut <- lapply(scope,function(var){
        varfit <- stats::update(object, as.formula(paste("~ . -", var)), 
                                evaluate = FALSE)
        varfit <- eval(varfit, envir = env)
        nnew <- stats::nobs(varfit, use.fallback = TRUE)
        if (all(is.finite(c(n0, nnew))) && nnew[[1]] != n0)
            stop("number of rows in use has changed: remove missing values?")
        varfit
    })
    names(leaveOneOut) <- scope
    if (missing(newdata)) newdata=eval(object$call$data)
    if (missing(formula)) formula=object$formula
    r2 <- Score(c("Full model"=list(object),leaveOneOut),
                formula=formula,
                data=newdata,
                contrasts=FALSE,
                conf.int=FALSE,
                metrics="brier",
                summary="rsquared",
                ...)$Brier$score
    ## r2[,IPA:=100*IPA]
    ## r2 <- r2[model!="Null model"]
    data.table::setnames(r2,"model","Variable")
    r2[,IPA.drop:=IPA[Variable=="Full model"]-IPA]
    ## r2 <- r2[Variable!="Full model"]
    out <- as.data.frame(r2[,data.table::data.table(Variable,Brier,IPA,IPA.drop)])
    class(out) <- c("IPA",class(out))
    out
} 

#' @export
IPA <- rsquared
#' @export
IPA.default <- rsquared.default
#' @export
IPA.CauseSpecificCox <- rsquared.CauseSpecificCox
#' @export
IPA.coxph <- rsquared.coxph
#' @export
IPA.glm <- rsquared.glm


######################################################################
### rsquared.R ends here
