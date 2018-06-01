### confint.influenceTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  1 2018 (12:15) 
## Version: 
## Last-Updated: jun  1 2018 (15:46) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.influenceTest (documentation)
##' @title Confidence Intervals and Confidence Bands for the Difference Between Two Estimates
##' @description Confidence intervals and confidence Bands for the difference between two estimates.
##' @name confint.influenceTest
##' 
##' @param object A \code{influenceTest} object, i.e. output of the \code{influenceTest} function.
##' @param parm not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] Level of confidence.
##' @param transform [character] the transformation used to improve coverage of the confidence intervals.
##' Can be \code{"none"} or \code{"atanh"}.
##' @param nsim.band [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set before performing simulations for the confidence bands.
##' If not given or NA no seed is set.
##' @param ... not used.
##'
##' @details Except for the cumulative hazard,
##' the confidence bands and confidence intervals are automatically restricted to the interval [-1;1].
##' 
##' @author Brice Ozenne
##'

## * confint.influenceTest (code)
##' @rdname confint.influenceTest
##' @export
confint.influenceTest <- function(object,
                                  parm = NULL,
                                  level = 0.95,
                                  nsim.band = 1e4,
                                  transform = "none",
                                  seed = NA,
                                  ...){

    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    if(!is.null(object$transform) && object$transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    if(object$type == "cumhazard"){
        object$transform <- match.arg(transform, c("none"))
    }else{
        object$transform <- match.arg(transform, c("none","atanh"))
    }

    ## ** compute se, CI/CB
    outCIBP <- transformCIBP(estimate = object$delta,
                             se = object$delta.se,
                             iid = object$delta.iid,
                             null = 0,
                             conf.level = level,
                             nsim.band = nsim.band,
                             seed = seed,
                             type = object$transform,
                             min.value = switch(object$transform,
                                                "none" = -1,
                                                "atanh" = NULL),
                             max.value = switch(object$transform,
                                                "none" = 1,
                                                "atanh" = 1),
                             ci = TRUE,
                             band = object$band,
                             p.value = TRUE)
    
    names(outCIBP) <- paste0("delta.", names(outCIBP))
    object[names(outCIBP)] <- outCIBP
    
    ## export
    object$conf.level <- level
    return(object)
}



######################################################################
### confint.influenceTest.R ends here
