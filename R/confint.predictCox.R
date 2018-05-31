### confint.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (18:05) 
##           By: Brice Ozenne
##     Update #: 198
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
##' @param nsim.band [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set when performing simulation for the confidence bands.
##' If not given or NA no seed is set.
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval of definition of the statistic,
##' i.e. a confidence interval for the survival of [0.5;1.2] will become [0.5;1].
##' 
##' 
##' @author Brice Ozenne
                                        
## * confint.predictCox (examples)
##' @rdname confint.predictCox
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
##'                        se = TRUE, band = TRUE)
##' fit.pred
##' sqrt(rowSums(fit.pred$survival.iid[1,,]^2)) ## se for individual 1
##'
##' #### add confidence intervals / bands computed on the original scale
##' confint(fit.pred, survival.transform = "none")
##' fit.pred$survival - 1.96 * fit.pred$survival.se  ## survival.lower
##' fit.pred$survival + 1.96 * fit.pred$survival.se  ## survival.upper
##'
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' confint(fit.pred, survival.transform = "loglog")
##' newse <- fit.pred$survival.se/(-fit.pred$survival*log(fit.pred$survival))
##' exp(-exp(log(-log(fit.pred$survival)) + 1.96 * newse)) ## survival.lower
##' exp(-exp(log(-log(fit.pred$survival)) - 1.96 * newse)) ## survival.upper
##'

## * confint.predictCox (code)
##' @rdname confint.predictCox
##' @export
confint.predictCox <- function(object,
                               level = 0.95,
                               parm = NULL,
                               nsim.band = 1e4,
                               cumhazard.transform = "log",
                               survival.transform = "cloglog",
                               seed = NA,
                               ...){

    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }
    
    if(is.null(parm)){
        parm <- intersect(c("survival","cumhazard"),names(object))[1]
    }
    if(length(parm)!=1){
        stop("Argument \'parm\' must have length 1 \n")
    }
    if(parm %in% c("cumhazard","survival") == FALSE){
        stop("Argument \'parm\' must be \"cumhazard\" or \"survival\" \n")
    }
    if(parm %in% names(object) == FALSE){
        stop(parm," has not been stored in the object \n",
             "set argument \'parm\' to ",parm," when calling predictCox \n")
    }
    if("cumhazard" %in% parm){
        if(!is.null(object$cumhazard.transform) && object$cumhazard.transform != "none"){
            stop("Cannot work with standard errors that have already been transformed \n")
        }
        object$cumhazard.transform <- match.arg(cumhazard.transform, c("none","log"))
        min.value <- switch(object$cumhazard.transform,
                            "none" = 0,
                            "log" = NULL)
        max.value <- NULL
    }
    if("survival" %in% parm){
        if(!is.null(object$survival.transform) && object$survival.transform != "none"){
            stop("Cannot work with standard errors that have already been transformed \n")
        }
        object$survival.transform <- match.arg(survival.transform, c("none","log","loglog","cloglog"))
        min.value <- switch(object$survival.transform,
                            "none" = 0,
                            "log" = NULL,
                            "loglog" = NULL,
                            "cloglog" = NULL)
        max.value <- switch(object$survival.transform,
                            "none" = 1,
                            "log" = 1,
                            "loglog" = NULL,
                            "cloglog" = NULL)
    }

    ## ** quantile
    zval <- stats::qnorm(1 - (1-level)/2, mean = 0, sd = 1)

    ## ** transformation
    if(object[[paste0(parm,".transform")]] != "none"){
        ## transform standard error
        object[[paste0(parm,".se")]] <- transformSE(estimate = object[[parm]],
                                                    se = object[[paste0(parm,".se")]],
                                                    type = object[[paste0(parm,".transform")]])

        ## transform influence function
        if(object$band){
            object[[paste0(parm,".iid")]] <- transformIID(estimate = object[[parm]],
                                                          iid = object[[paste0(parm,".iid")]],
                                                          type = object[[paste0(parm,".transform")]],
                                                          format = "array")
        }
    }

    ## ** confidence intervals
    if(object$se){

        object[paste0(parm,c(".lower",".upper"))] <- transformCI(estimate = object[[parm]],
                                                                 se = object[[paste0(parm,".se")]],
                                                                 quantile = zval,
                                                                 type = object[[paste0(parm,".transform")]],
                                                                 format = "matrix",
                                                                 min.value = min.value,
                                                                 max.value = max.value)

    }

    ## ** confidence bands
    if(object$band && nsim.band > 0){

        ## find quantiles for the bands
        if(!is.na(seed)){set.seed(seed)}
        ## sqrt(rowSums(object[[paste0(parm,".iid")]][1,,]^2))
        object[[paste0(parm,".quantileBand")]] <- confBandCox(iid = object[[paste0(parm,".iid")]],
                                                              se = object[[paste0(parm,".se")]],
                                                              n.sim = nsim.band,
                                                              conf.level = level)

        object[paste0(parm,c(".lowerBand",".upperBand"))] <- transformCI(estimate = object[[parm]],
                                                                         se = object[[paste0(parm,".se")]],
                                                                         quantile = object[[paste0(parm,".quantileBand")]],
                                                                         type = object[[paste0(parm,".transform")]],
                                                                         format = "matrix",
                                                                         min.value = min.value,
                                                                         max.value = max.value)

    }

    ## ** check NA
    indexNA <- union(which(is.na(object[[paste0(parm,".se")]])),
                     which(is.nan(object[[paste0(parm,".se")]])))
    if(length(indexNA)>0){

        if(object$se){
            object[[paste0(parm,".lower")]][indexNA] <- NA
            object[[paste0(parm,".upper")]][indexNA] <- NA
        }
        if(object$band){ ## if cannot compute se at one time then remove confidence band at all times
            indexNA2 <- union(which(rowSums(is.na(object[[paste0(parm,".se")]]))>0),
                              which(rowSums(is.nan(object[[paste0(parm,".se")]]))>0))
            object[[paste0(parm,".quantileBand")]][indexNA2] <- NA
            object[[paste0(parm,".lowerBand")]][indexNA2,] <- NA
            object[[paste0(parm,".upperBand")]][indexNA2,] <- NA
        }
        
    }

    
    ## ** export
    object$conf.level <- level
    return(object)
}

######################################################################
### confint.predictCox.R ends here
