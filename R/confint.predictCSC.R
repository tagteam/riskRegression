### confint.predictCSC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (10:31) 
##           By: Brice Ozenne
##     Update #: 123
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.predictCSC (documentation)
##' @title Confidence Intervals and Confidence Bands for the Predicted Absolute Risk (Cumulative Incidence Function)
##' @description Confidence intervals and confidence Bands for the predicted absolute risk (cumulative incidence function).
##' @name confint.predictCSC
##' 
##' @param object A \code{predictCSC} object, i.e. output of the \code{predictCSC} function.
##' @param conf.level Level of confidence.
##' @param absRisk.transform the transformation used to improve coverage
##' of the confidence intervals for the predicted absolute risk in small samples.
##' @param nsim.band the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed Integer passed to set.seed when performing simulation for the confidence bands.
##' If not given or NA no seed is set. 
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval [0;1].
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

## * confint.predictCSC (code)
##' @rdname confint.predictCSC
##' @export
confint.predictCSC <- function(object,
                               conf.level = 0.95,
                               nsim.band = 1e4,
                               absRisk.transform = "cloglog",
                               seed = NA,
                               ...){

    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    if(!is.null(object$absRisk.transform) && object$absRisk.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    object$absRisk.transform <- match.arg(absRisk.transform, c("none","log","loglog","cloglog"))

    ## ** quantile
    zval <- stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1)

    ## ** transformation
    if(object$absRisk.transform != "none"){
        ## transform standard error
        object$absRisk.se <- transformSE(estimate = object$absRisk,
                                         se = object$absRisk.se,
                                         type = object$absRisk.transform)

        ## transform influence function
        if(object$band){
            object$absRisk.iid <- transformIID(estimate = object$absRisk,
                                               iid = object$absRisk.iid,
                                               type = object$absRisk.transform,
                                               format = "array")
        }

        min.value <- NULL
        max.value <- NULL
    }

    min.value <- switch(object$absRisk.transform,
                        "none" = 0,
                        "log" = NULL,
                        "loglog" = NULL,
                        "cloglog" = NULL)
    max.value <- switch(object$absRisk.transform,
                        "none" = 1,
                        "log" = 1,
                        "loglog" = NULL,
                        "cloglog" = NULL)
    
    ## ** confidence intervals
    if(object$se){

        object[c("absRisk.lower","absRisk.upper")] <- transformCI(estimate = object$absRisk,
                                                                  se = object$absRisk.se,
                                                                  quantile = zval,
                                                                  type = object$absRisk.transform,
                                                                  format = "matrix",
                                                                  min.value = min.value,
                                                                  max.value = max.value)
    }
        

    ## confidence bands
    if(object$band && nsim.band > 0){

        ## find quantiles for the bands
        if(!is.na(seed)){set.seed(seed)}
        object$absRisk.quantileBand <- confBandCox(iid = object$absRisk.iid,
                                                   se = object$absRisk.se,
                                                   n.sim = nsim.band,
                                                   conf.level = conf.level)

        object[c("absRisk.lowerBand","absRisk.upperBand")] <- transformCI(estimate = object$absRisk,
                                                                          se = object$absRisk.se,
                                                                          quantile = object$absRisk.quantileBand,
                                                                          type = object$absRisk.transform,
                                                                          format = "matrix",
                                                                          min.value = min.value,
                                                                          max.value = max.value)
    }

    ## check NA
    indexNA <- union(which(is.na(object$absRisk.se)),
                     which(is.nan(object$absRisk.se)))
    if(length(indexNA)>0){

        if(se){
            object$absRisk.lower[indexNA] <- NA
            object$absRisk.upper[indexNA] <- NA
        }
        if(object$band){
            indexNA2 <- union(which(rowSums(is.na(object$absRisk.se))>0),
                              which(rowSums(is.nan(object$absRisk.se))>0))
            object$absRisk.quantileBand[indexNA2] <- NA
            object$absRisk.lowerBand[indexNA2,] <- NA
            object$absRisk.upperBand[indexNA2,] <- NA
        }
        
    }

    
    ## export
    object$conf.level <- conf.level
    return(object)
}


######################################################################
### confint.predictCSC.R ends here
