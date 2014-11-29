# inverse of the probability of censoring weigths at the subject specific event times
# {{{ method root
#' Estimation of censoring probabilities at subject specific times
#' 
#' This function is used internally to contruct pseudo values by inverse of the
#' probability of censoring weights.
#' 
#' Inverse of the probability of censoring weights usually refer to the
#' probabilities of not being censored at certain time points. These
#' probabilities are also the values of the conditional survival function of
#' the censoring time given covariates. The function subjectWeights estimates
#' the conditional survival function of the censoring times and derives the
#' weights.
#' 
#' IMPORTANT: the data set should be ordered, \code{order(time,-status)} in
#' order to get the \code{weights} in the right order for some choices of
#' \code{method}.
#' 
#' @aliases subjectWeights subjectWeights.none subjectWeights.km
#' subjectWeights.marginal subjectWeights.nonpar subjectWeights.cox
#' subjectWeights.aalen
#' @param formula A survival formula like, Surv(time,status)~1 or
#' Hist(time,status)~1 where status=0 means censored. The status variable is
#' internally reversed for estimation of censoring rather than survival
#' probabilities. Some of the available models, see argument \code{model}, will
#' use predictors on the right hand side of the formula.
#' @param data The data used for fitting the censoring model
#' @param method Censoring model used for estimation of the (conditional)
#' censoring distribution.
#' @param lag If equal to \code{1} then obtain \code{G(T_i-|X_i)}, if equal to
#' \code{0} estimate the conditional censoring distribution at the
#' subjectTimes, i.e. (\code{G(T_i|X_i)}).
#' @return \item{times}{The times at which weights are estimated}
#' \item{weights}{Estimated weights at individual time values
#' \code{subjectTimes}} \item{lag}{The time lag.} \item{fit}{The fitted
#' censoring model} \item{method}{The method for modelling the censoring
#' distribution} \item{call}{The call}
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @keywords survival
#' @examples
#' 
#' library(pec)
#' dat=SimSurv(300)
#' 
#' dat <- dat[order(dat$time,-dat$status),]
#' 
#' # using the marginal Kaplan-Meier for the censoring times
#' 
#' WKM=subjectWeights(Hist(time,status)~X2,data=dat,method="marginal")
#' plot(WKM$fit)
#' WKM$fit
#' WKM$weights
#' 
#' # using the Cox model for the censoring times given X2
#' 
#' WCox=subjectWeights(Surv(time,status)~X2,data=dat,method="cox")
#' WCox
#' plot(WCox$weights,WKM$weights)
#' 
#' # using the stratified Kaplan-Meier for the censoring times given X2
#' 
#' WKM2 <- subjectWeights(Surv(time,status)~X2,data=dat,method="nonpar")
#' plot(WKM2$fit,add=FALSE)
#' 
#'
#' @export subjectWeights
subjectWeights <- function(formula,data,method=c("cox","marginal","km","nonpar","rfsrc","none"),args,lag=1){
    if (lag!=1 && lag!=0){ stop("lag must be either 0 or 1")}
    method <- tolower(method)
    method <- match.arg(method,c("cox","marginal","km","nonpar","rfsrc","none"))
    class(method) <- method
    UseMethod("subjectWeights",method)
}
# }}}
# {{{ None: set weights to 1
#' @S3method subjectWeights none
subjectWeights.none <- function(formula,data,method,args,lag=1){
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
#' @S3method subjectWeights marginal
subjectWeights.marginal <- function(formula,data,method,args,lag=1){
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
#' @S3method subjectWeights nonpar
subjectWeights.nonpar <- function(formula,data,method,args,lag=1){
    fit <- prodlim(formula,data=data,reverse=TRUE,bandwidth="smooth")
    weights <- predictSurvIndividual(fit,lag=lag)
    out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
    class(out) <- "subjectWeights"
    out
}
# }}}
# {{{ reverse Cox via Harrel's package
#' @S3method subjectWeights cox
subjectWeights.cox <- function(formula,data,method,args,lag=1){
    ## require(rms)
    require(survival)
    EHF <- prodlim::EventHistory.frame(formula,data,specials=NULL,unspecialsDesign=FALSE)
    wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    stopifnot(NROW(na.omit(wdata))>0)    
    fit <- cph(formula,data=wdata,surv=TRUE,x=TRUE,y=TRUE)
    if (lag==1)
        weights <- as.vector(survest(fit,times=times-min(diff(c(0,unique(times))))/2,what='parallel'))
    else # (lag==0)
        weights <- as.vector(survest(fit,times=times,what='parallel'))
    out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
    class(out) <- "subjectWeights"
    out
}
# }}}
# {{{ reverse random forest
#' @S3method subjectWeights rfsrc
subjectWeights.rfsrc <- function(formula,data,method,args,lag=1){
    EHF <- prodlim::EventHistory.frame(formula,data,specials=NULL,unspecialsDesign=FALSE)
    wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    stopifnot(NROW(na.omit(wdata))>0)    
    if (missing(args) || is.null(args))
        ## args <- list(bootstrap="none",ntree=1000,nodesize=NROW(data)/2)
        args <- list(ntree=1000)
    ## if (is.null(args$importance) & (args$importance!="none"))
    args$importance <- "none"
    fit <- do.call(randomForestSRC::rfsrc,c(list(wform,data=wdata),args))
    pmat <- fit$survival
    jtimes <- fit$time.interest
    IPCW.subjectTimes <- sapply(1:length(subjectTimes),function(i){
        Ci <- subjectTimes[i]
        pos <- prodlim::sindex(jump.times=jtimes,eval.times=Ci,comp="smaller",strict=(lag==1))
        c(1,pmat[i,])[1+pos]
    })
    weights <- as.vector(survest(fit,times=times,what='parallel'))
    out <- list(weights=weights,fit=fit,lag=lag,call=match.call(),method=method)
    class(out) <- "subjectWeights"
    out
}
# }}}
