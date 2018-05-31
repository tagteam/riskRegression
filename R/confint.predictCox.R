### confint.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (09:40) 
##           By: Brice Ozenne
##     Update #: 177
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
##' @param conf.level Level of confidence.
##' @param type the type of predicted value for which the confidence intervals should be output.
##' Can be \code{"survival"} or \code{"cumhazard"}.
##' @param cumhazard.transform the transformation used to improve coverage
##' of the confidence intervals for the cumlative hazard in small samples.
##' @param survival.transform the transformation used to improve coverage
##' of the confidence intervals for the survival in small samples.
##' @param nsim.band the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed Integer passed to set.seed when performing simulation for the confidence bands.
##' If not given or NA no seed is set. 
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval of definition of the statistic,
##' i.e. a confidence interval for the survival of [0.5;1.2] will become [0.5;1].
##' 
##' 
##' @author Brice Ozenne
##'
##' @examples
##' library(survival)
##'
##' ## generate data
##' set.seed(10)
##' d <- sampleData(40,outcome="survival") ## training dataset
##' nd <- sampleData(4,outcome="survival") ## validation dataset
##' d$time <- round(d$time,1) ## create tied events
##' # table(duplicated(d$time))
##' 
##' ## estimate a stratified Cox model
##' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
##'              data=d, ties="breslow", x = TRUE, y = TRUE)
##'
##' 
##' ## compute individual specific cumulative hazard and survival probabilities 
##' fit.pred <- predictCox(fit, newdata=nd, times=c(3,8), se = TRUE)
##' fit.pred
##'
##' ## add confidence intervals computed on the original scale
##' confint(fit.pred, transform.survival = "none")
##' fit.pred$survival + 1.96 * fit.pred$survival.se  ## survival.upper
##'
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' confint(fit.pred, transform.survival = "loglog")
##' newse <- fit.pred$survival.se/(fit.pred$survival*log(fit.pred$survival))
##' exp(-exp(log(-log(fit.pred$survival)) - 1.96 * newse)) ## survival.lower
##' exp(-exp(log(-log(fit.pred$survival)) + 1.96 * newse)) ## survival.upper
##

## * confint.predictCox (code)
##' @rdname confint.predictCox
##' @export
confint.predictCox <- function(object,
                               conf.level = 0.95,
                               type = NULL,
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
    
    if(is.null(type)){
        type <- intersect(c("survival","cumhazard"),names(object))[1]
    }
    if(length(type)!=1){
        stop("Argument \'type\' must have length 1 \n")
    }
    if(type %in% c("cumhazard","survival") == FALSE){
        stop("Argument \'type\' must be \"cumhazard\" or \"survival\" \n")
    }
    if(type %in% names(object) == FALSE){
        stop(type," has not been stored in the object \n",
             "set argument \'type\' to ",type," when calling predictCox \n")
    }
    if("cumhazard" %in% type){
        if(!is.null(object$cumhazard.transform) && object$cumhazard.transform != "none"){
            stop("Cannot work with standard errors that have already been transformed \n")
        }
        object$cumhazard.transform <- match.arg(cumhazard.transform, c("none","log"))
        min.value <- switch(object$cumhazard.transform,
                            "none" = 0,
                            "log" = NULL)
        max.value <- NULL
    }
    if("survival" %in% type){
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
    zval <- stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1)

    ## ** transformation
    if(object[[paste0(type,".transform")]] != "none"){
        ## transform standard error
        object[[paste0(type,".se")]] <- transformSE(estimate = object[[type]],
                                                    se = object[[paste0(type,".se")]],
                                                    type = object[[paste0(type,".transform")]])

        ## transform influence function
        if(object$band){
            object[[paste0(type,".iid")]] <- transformIID(estimate = object[[type]],
                                                          iid = object[[paste0(type,".iid")]],
                                                          type = object[[paste0(type,".transform")]],
                                                          format = "array")
        }
    }
    
    ## ** confidence intervals
    if(object$se){

        object[paste0(type,c(".lower",".upper"))] <- transformCI(estimate = object[[type]],
                                                                 se = object[[paste0(type,".se")]],
                                                                 quantile = zval,
                                                                 type = object[[paste0(type,".transform")]],
                                                                 format = "matrix",
                                                                 min.value = min.value,
                                                                 max.value = max.value)

    }

    ## ** confidence bands
    if(object$band && nsim.band > 0){

        ## find quantiles for the bands
        if(!is.na(seed)){set.seed(seed)}        
        object[[paste0(type,".quantileBand")]] <- confBandCox(iid = object[[paste0(type,".iid")]],
                                                              se = object[[paste0(type,".se")]],
                                                              n.sim = nsim.band,
                                                              conf.level = conf.level)

        object[paste0(type,c(".lowerBand",".upperBand"))] <- transformCI(estimate = object[[type]],
                                                                         se = object[[paste0(type,".se")]],
                                                                         quantile = object[[paste0(type,".quantileBand")]],
                                                                         type = object[[paste0(type,".transform")]],
                                                                         format = "matrix",
                                                                         min.value = min.value,
                                                                         max.value = max.value)

    }

    ## ** check NA
    indexNA <- union(which(is.na(object[[paste0(type,".se")]])),
                     which(is.nan(object[[paste0(type,".se")]])))
    if(length(indexNA)>0){

        if(object$se){
            object[[paste0(type,".lower")]][indexNA] <- NA
            object[[paste0(type,".upper")]][indexNA] <- NA
        }
        if(object$band){ ## if cannot compute se at one time then remove confidence band at all times
            indexNA2 <- union(which(rowSums(is.na(object[[paste0(type,".se")]]))>0),
                              which(rowSums(is.nan(object[[paste0(type,".se")]]))>0))
            object[[paste0(type,".quantileBand")]][indexNA2] <- NA
            object[[paste0(type,".lowerBand")]][indexNA2,] <- NA
            object[[paste0(type,".upperBand")]][indexNA2,] <- NA
        }
        
    }

    
    ## ** export
    object$conf.level <- conf.level
    return(object)
}

######################################################################
### confint.predictCox.R ends here
