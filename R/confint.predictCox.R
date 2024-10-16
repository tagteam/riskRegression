### confint.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: Oct 15 2024 (11:48) 
##           By: Brice Ozenne
##     Update #: 345
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.predictCox (documentation)
##' @title Confidence Intervals and Confidence Bands for the predicted Survival/Cumulative Hazard
##' @description Confidence intervals and confidence Bands for the predicted survival/cumulative Hazard.
##' @name confint.predictCox
##' 
##' @param object A \code{predictCox} object, i.e. output of the \code{predictCox} function.
##' @param level [numeric, 0-1] Level of confidence.
##' @param parm [character] the type of predicted value for which the confidence intervals should be output.
##' Can be \code{"survival"} or \code{"cumhazard"}.
##' @param cumhazard.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the cumlative hazard in small samples.
##' Can be \code{"none"}, \code{"log"}.
##' @param survival.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the survival in small samples.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}.
##' @param n.sim [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set before performing simulations for the confidence bands.
##' If not given or NA no seed is set.
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval of definition of the statistic,
##' i.e. a confidence interval for the survival of [0.5;1.2] will become [0.5;1].
##' 
##' 
##' @author Brice Ozenne
                                        
## * confint.predictCox (examples)
##' @examples
##' library(survival)
##'
##' #### generate data ####
##' set.seed(10)
##' d <- sampleData(40,outcome="survival") 
##' 
##' #### estimate a stratified Cox model ####
##' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
##'              data=d, ties="breslow", x = TRUE, y = TRUE)
##' 
##' #### compute individual specific survival probabilities  
##' fit.pred <- predictCox(fit, newdata=d[1:3], times=c(3,8), type = "survival",
##'                        se = TRUE, iid = TRUE, band = TRUE)
##' fit.pred
##'
##' ## check standard error
##' sqrt(rowSums(fit.pred$survival.iid[,,1]^2)) ## se for individual 1
##'
##' ## check confidence interval
##' newse <- fit.pred$survival.se/(-fit.pred$survival*log(fit.pred$survival))
##' cbind(lower = as.double(exp(-exp(log(-log(fit.pred$survival)) + 1.96 * newse))),
##'       upper = as.double(exp(-exp(log(-log(fit.pred$survival)) - 1.96 * newse)))
##' )
##'
##' #### compute confidence intervals without transformation
##' confint(fit.pred, survival.transform = "none")
##' cbind(lower = as.double(fit.pred$survival - 1.96 * fit.pred$survival.se),
##'       upper = as.double(fit.pred$survival + 1.96 * fit.pred$survival.se)
##' )
##'

## * confint.predictCox (code)
##' @rdname confint.predictCox
##' @method confint predictCox
##' @export
confint.predictCox <- function(object,
                               parm = NULL,
                               level = 0.95,
                               n.sim = 1e4,
                               cumhazard.transform = "log",
                               survival.transform = "loglog",
                               seed = NA,
                               ...){

    if(object$se[[1]] == FALSE && object$band[[1]] == FALSE){
        message("No confidence interval/band computed \n",
                "Set argument \'se\' or argument \'band\' to TRUE when calling the predictCox function \n")
        return(object)
    }

    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }
    
    if(is.null(parm)){
        parm <- intersect(c("lp","survival","cumhazard"),names(object))
    }else if(any(parm %in% c("cumhazard","survival") == FALSE)){
        txt <- parm[parm %in% c("cumhazard","survival") == FALSE]
        txt2 <- paste0("\"",paste0(txt, collapse = "\" \""),"\"")
        stop("Argument \'parm\' must be \"cumhazard\" or \"survival\" \n",
             "incorrect value(s): ",txt2," \n")
    }

    if(any(parm %in% names(object) == FALSE)){
        txt <- parm[parm %in% names(object) == FALSE]
        txt2 <- paste0("\"",paste0(txt, collapse = "\" \""),"\"")
        stop(txt2," has/have not been stored in the object \n",
             "set argument \'parm\' to ",txt2," when calling the predictCox function \n")
    }

    if("lp" %in% parm){
        object$lp.transform <- "none"
    }
    
    if("cumhazard" %in% parm){
        object$cumhazard.transform <- match.arg(cumhazard.transform, c("none","log"))
        if(object$band[[1]] && is.null(object$cumhazard.se)){
            stop("Cannot compute confidence bands \n",
                 "Set argument \'se\' to TRUE when calling the predictCox function \n")
        }
        if(object$band[[1]] && is.null(object$cumhazard.iid)){
            stop("Cannot compute confidence bands \n",
                 "Set argument \'iid\' to TRUE when calling the predictCox function \n")
        }

    }
    
    if("survival" %in% parm){
        object$survival.transform <- match.arg(survival.transform, c("none","log","loglog","cloglog"))
        if(object$band[[1]] && is.null(object$survival.se)){
            stop("Cannot compute confidence bands \n",
                 "Set argument \'se\' to TRUE when calling the predictCox function \n")
        }
        if(object$band[[1]] && is.null(object$survival.iid)){
            stop("Cannot compute confidence bands \n",
                 "Set argument \'iid\' to TRUE when calling the predictCox function \n")
        }
    }

    ## ** compute se, CI/CB
    object$vcov <- setNames(vector(mode = "list", length = length(parm)), parm)
    for(iType in parm){
        if(iType=="lp"){
            iMin.value <- NULL
            iMax.value <- NULL
            iEstimate <- matrix(object[[iType]], nrow = 1)
            if(object$se[[1]]){
                iSe <- matrix(object[[paste0(iType,".se")]], nrow = 1)
            }else{
                iSe <- NULL
            }
            if(object$band[[1]]){
                iIID <- array(object[[paste0(iType,".iid")]], dim = c(NROW(object[[paste0(iType,".iid")]]),NCOL(object[[paste0(iType,".iid")]]),1))
            }else{
                iIID <- NULL
            }
        }else if(iType=="cumhazard"){
            iMin.value <- switch(object$cumhazard.transform,
                                 "none" = 0,
                                 "log" = NULL)
            iMax.value <- NULL
            iEstimate <- object[[iType]]
            iSe <- object[[paste0(iType,".se")]]
            iIID <- object[[paste0(iType,".iid")]]
        }else if(iType=="survival"){
            iMin.value <- switch(object$survival.transform,
                                 "none" = 0,
                                 "log" = NULL,
                                 "loglog" = NULL,
                                 "cloglog" = NULL)
            iMax.value <- switch(object$survival.transform,
                                 "none" = 1,
                                 "log" = 1,
                                 "loglog" = NULL,
                                 "cloglog" = NULL)
            iEstimate <- object[[iType]]
            iSe <- object[[paste0(iType,".se")]]
            iIID <- object[[paste0(iType,".iid")]]
        }

        ## reshape to multiple adjust across subject instead of timepoint when using diag
        if(object$diag && object$band){
            iEstimate <- t(iEstimate)
            iSe <- t(iSe)
            iIID <- aperm(iIID, c(1,3,2))
        }
        outCIBP <- transformCIBP(estimate = iEstimate,
                                 se = iSe,
                                 iid = iIID,
                                 null = NA,
                                 conf.level = level,
                                 n.sim = n.sim,
                                 seed = seed,
                                 type = object[[paste0(iType,".transform")]],
                                 min.value = iMin.value,
                                 max.value = iMax.value,
                                 ci = object$se,
                                 band = object$band,
                                 method.band = "maxT-simulation",
                                 alternative = "two.sided",
                                 p.value = FALSE)
        ## restaure original shape
        if(object$diag && object$band){
            outCIBP$lower <- t(outCIBP$lower)
            outCIBP$upper <- t(outCIBP$upper)
            outCIBP$lowerBand <- t(outCIBP$lowerBand)
            outCIBP$upperBand <- t(outCIBP$upperBand)
        }
        names(outCIBP) <- paste0(iType,".", names(outCIBP))
        object[names(outCIBP)] <- outCIBP
        
        ## ** restaure dimensions
        if(iType=="lp"){
            if(object$se[[1]]){
                object$lp.lower <- matrix(object$lp.lower, ncol = 1)
                object$lp.upper <- matrix(object$lp.upper, ncol = 1)
            }
            if(object$band[[1]]){
                object$lp.lowerBand <- matrix(object$lp.lowerBand, ncol = 1)
                object$lp.upperBand <- matrix(object$lp.upperBand, ncol = 1)
            }
        }
    }

    ## ** export
    object$conf.level <- level
    return(object)
}

######################################################################
### confint.predictCox.R ends here
