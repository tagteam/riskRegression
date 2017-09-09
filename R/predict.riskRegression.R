#' Predict individual risk.
#' 
#' Extract predictions from a risk prediction model.
#' 
#' 
#' @param object Fitted object obtained with one of \code{ARR}, \code{LRR},
#' \code{riskRegression}.
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted risk.
#' @param \dots not used
#' @author Thomas H. Scheike \email{ts@@biostat.ku.dk}
#' 
#' Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @references Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute risk
#' regression for competing risks: interpretation, link functions and
#' prediction Research report 11/8. Department of Biostatistics, University of
#' Copenhagen
#' @keywords survival
##' @examples
##' 
##' data(Melanoma)
##' library(prodlim)
##' library(survival)
##' 
##' fit.tarr <- ARR(Hist(time,status)~age+invasion+strata(sex),data=Melanoma,cause=1)
##' predict(fit.tarr,newdata=data.frame(age=48,
##'                      invasion=factor("level.1",
##'                          levels=levels(Melanoma$invasion)),
##'                      sex=factor("Female",levels=levels(Melanoma$sex))))
##' predict(fit.tarr,newdata=data.frame(age=48,
##'                      invasion=factor("level.1",
##'                          levels=levels(Melanoma$invasion)),
##'                      sex=factor("Male",levels=levels(Melanoma$sex))))
##' predict(fit.tarr,newdata=data.frame(age=c(48,58,68),
##'                      invasion=factor("level.1",
##'                          levels=levels(Melanoma$invasion)),
##'                      sex=factor("Male",levels=levels(Melanoma$sex))))
##' predict(fit.tarr,newdata=Melanoma[1:4,])
##'
#' @method predict riskRegression
#' @export 
predict.riskRegression <- function(object,
                                   newdata,
                                   ...){
    Zcoef <- c(object$timeConstantEffects$coef)
    semi <- !is.null(Zcoef)
    tt <- delete.response(object$design$Terms)
    PF <- prodlim::EventHistory.frame(formula(tt),
                                      newdata,
                                      specials=c("timevar","strata","prop","const","tp"),
                                      stripSpecials=c("timevar","prop"),
                                      stripArguments=list("prop"=list("power"=0),"timevar"=list("test"=0)),
                                      stripAlias=list("timevar"=c("strata"),"prop"=c("tp","const")),
                                      stripUnspecials="prop",
                                      specialsDesign=TRUE,
                                      dropIntercept=TRUE,check.formula=FALSE,response=FALSE)
    ##  The time-constant effects
    Z <- PF$prop
    factorLevelsZ <- attr(PF$prop,"levels")
    refLevelsZ <- lapply(factorLevelsZ,function(x)x[1])
    if (!(all(colnames(Z) %in% names(object$design$timepower))))
        stop("\nProblem with factor names.\nCheck if the storage type of all factors\nin the original data are compatible with those in newdata\n\nOffending variable(s): ",paste(colnames(Z)[!colnames(Z)%in% names(object$timePower)],collapse=", ")," does not match ",paste(names(object$timePower)[!names(object$timePower)%in%colnames(Z)],collapse=" ,"),".")
    ## The time-varying effects
    X <- PF$timevar
    factorLevelsX <- attr(PF$timevar,"levels")
    refLevelsX <- lapply(factorLevelsX,function(x)x[1])
    ## remove the time column 
    Xcoef  <-  object$timeVaryingEffects$coef[,-1,drop=FALSE]
    # {{{ compute the linear predictor
    fittime  <-  object$time
    ntime  <-  NROW(fittime)
    ## intercept
    if (NROW(X)==0){
        ## X <- cbind("Intercept"=1)
        ## remove the time column
        intercept <- t(Xcoef)
        timeVarLP <- matrix(rep(intercept,NROW(Z)),
                            byrow=TRUE,
                            nrow=NROW(Z))
    }
    else {
        X <- cbind("Intercept"=1,X)
        colnamesX <- c("Intercept",colnames(X))
        timeVarLP <- X %*% t(Xcoef) 
        ## FIXME: manually set extreme values if LP = 0
        ##        should be done in the c-routine
        timeVarLP <- t(apply(timeVarLP,1,function(lp){
            fixedlp <- lp;
            fixedlp[fixedlp==0] <- -Inf
            fixedlp}))
    }
    if (semi==TRUE) {
        timepower <- object$design$timepower
        if (any(timepower>0)){
            timeFactor <- matrix(fittime,nrow=length(Zcoef),ncol=length(fittime),byrow=TRUE)
            timeConstLP  <-  Z %*% (Zcoef * (timeFactor^timepower))
        }
        else{
            timeConstLP  <-  colSums(t(Z) * Zcoef)
        }
    }
    else
        timeConstLP <- 0
    LP <- timeConstLP+timeVarLP
    ## tag
    ## added 23 Jan 2016 (16:33)
    ## to fix bug in predictRisk.riskRegression
    if (length(fittime)==1)
        LP <- t(LP)
    # }}}
    # {{{ compute P1
    P1 <- switch(object$link,"relative"={
        pmin(exp(LP),1)
    },"additive"={
        P1=1-exp(-LP)
    },"prop"={
        P1 <- 1-exp(-exp(LP))
    },"logistic"={
        P1 <- exp(LP)/(1+exp(LP))
    })
    # }}}
    # {{{ standard errors and confidence bands
    ## use timereg
    # }}}
    # {{{ output
    out <- list(time=fittime,risk=P1)
    class(out) <- "predictedRisk"
    out
    # }}}
}
