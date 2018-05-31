### confint.predictCSC.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (16:58) 
##           By: Brice Ozenne
##     Update #: 130
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
##' @param parm not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] Level of confidence.
##' @param absRisk.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the predicted absolute risk in small samples.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}.
##' @param nsim.band [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set when performing simulation for the confidence bands.
##' If not given or NA no seed is set.
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval [0;1].
##' 
##' @author Brice Ozenne
##'

## * confint.predictCSC (examples)
##' @rdname confint.predictCSC
##' @examples
##' library(survival)
##'
##' #### generate data ####
##' set.seed(10)
##' d <- sampleData(100) 
##' 
##' #### estimate a stratified CSC model ###
##' fit <- CSC(Hist(time,event)~ X1 + strata(X2) + X6, data=d)
##' 
##' #### compute individual specific risks
##' fit.pred <- predict(fit, newdata=d[1:3], times=c(3,8), cause = 1,
##'                     se = TRUE, band = TRUE)
##' fit.pred
##'
##' ## add confidence intervals computed on the original scale
##' confint(fit.pred, absRisk.transform = "none")
##' fit.pred$absRisk - 1.96 * fit.pred$absRisk.se  ## survival.lower
##' fit.pred$absRisk + 1.96 * fit.pred$absRisk.se  ## survival.upper
##'
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' confint(fit.pred, absRisk.transform = "loglog")
##' newse <- fit.pred$absRisk.se/(-fit.pred$absRisk*log(fit.pred$absRisk))
##' exp(-exp(log(-log(fit.pred$absRisk)) + 1.96 * newse)) ## survival.lower
##' exp(-exp(log(-log(fit.pred$absRisk)) - 1.96 * newse)) ## survival.upper

## * confint.predictCSC (code)
##' @rdname confint.predictCSC
##' @export
confint.predictCSC <- function(object,
                               parm = NULL,
                               level = 0.95,
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
    zval <- stats::qnorm(1 - (1-level)/2, mean = 0, sd = 1)

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
                                                   conf.level = level)

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

        if(object$se){
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
    object$conf.level <- level
    return(object)
}


######################################################################
### confint.predictCSC.R ends here
