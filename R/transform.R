### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 30 2018 (15:58) 
## Version: 
## Last-Updated: maj 31 2018 (18:19) 
##           By: Brice Ozenne
##     Update #: 79
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * transformSE
##' @title Compute the Standard Error after Transformation
##' @description Compute the standard error after transformation
##' based on the standard error before transformation.
##'
##' @param estimate [numeric vector/matrix] the estimate value before transformation.
##' @param se [numeric vector/matrix] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' 
##' @details Use a delta method to find the standard error after transformation.
##'
##' @export
transformSE <- function(estimate, se, type){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.

    if(type == "none"){
        ## no change
        return(se)
    }else if(type == "log"){
        ## formula 4.10 p 58 (Beyersmann et al. 2012)
        return(se/estimate)
    }else if(type == "loglog"){
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        return(se / ( estimate * (-log(estimate)) ))
    }else if(type == "cloglog"){
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
        return(se / ( (1-estimate) * (-log(1-estimate)) ))
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
        ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
        ##                         = dx/(1+x^2)
        ##               Var(f(x)) = Var(x)/(1+x^2)
        return(se / (1+estimate^2) )
    }
}

## * transformIID
##' @title Compute the Influence Function after Transformation
##' @description Compute the influence function after transformation
##' based on the influence function before transformation.
##'
##' @param estimate [numeric vector/matrix] the estimate value before transformation.
##' @param iid [numeric matrix/array] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @param format [character] the class of iid, either \code{"matrix"} or \code{"array"}.
##' 
##' @details Use a delta method to find the standard error after transformation. \cr \cr
##'
##' Matrix case: the iid decomposition must contain have dimension [n.prediction,n.obs] and estimate [n.prediction]. \cr
##' Array case: the iid decomposition must contain have dimension [n.prediction,time,n.obs] and estimate [n.prediction,time].
##'
##' @export
transformIID <- function(estimate, iid, type, format){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.

    if(type == "none"){
        ## no change
        return(iid)
    }else if(type == "log"){
        ## formula 4.10 p 58 (Beyersmann et al. 2012)
        scaling <- estimate
    }else if(type == "loglog"){
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        scaling <- estimate * (-log(estimate))
    }else if(type == "cloglog"){
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
        scaling <- (1-estimate) * (-log(1-estimate))
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
        ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
        ##                         = dx/(1+x^2)
        ##               Var(f(x)) = Var(x)/(1+x^2)
        scaling <- 1 / (1+estimate^2)
    }

    if(format == "matrix"){
        return(colScale_cpp(iid, scale = scaling))
    }else if(format == "array"){
        return(sliceScale_cpp(iid, M = scaling))
    }
}

## * transformCI
##' @title Compute the Confidence Interval using a transformation
##' @description Compute the confidence interval using a transformation.
##' The resulting confidence interval is returned on the original case (i.e. back-transformed).
##'
##' @param estimate [numeric vector/matrix] the estimate value before transformation.
##' @param se [numeric matrix/array] the standard error after transformation.
##' @param quantile [numeric] quantile that will be multiplied to the standard error.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @param format [character] the class of iid, either \code{"matrix"} or \code{"array"}.
##' @param min.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{min},
##' it will be set at \code{min}. 
##' @param max.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{max},
##' it will be set at \code{max}.
##'
##' @details Argument \code{quantile} can either be a numeric value or a vector of numeric values.
##' In the latter case, if \code{se} is a matrix, \code{se} will be multiplied by \code{quantile} by column.
##' @export
transformCI <- function(estimate, se, quantile, type, format, min.value, max.value){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.
    out <- list()

    if(format == "vector" || length(quantile) == 1){
        quantileSe <- quantile * se
    }else{
        quantileSe <- colMultiply_cpp(se, scale = quantile)
    }
    
    ## compute confidence intervals
    if(type == "none"){
        out$lower <- estimate - quantileSe
        out$upper <- estimate + quantileSe
    }else if(type == "log"){
        ## a * exp +/-b = exp(log(a) +/- b)
        ## formula 4.10 p 58 (Beyersmann et al. 2012)
        out$lower <- estimate * exp(- quantileSe)
        out$upper <- estimate * exp(+ quantileSe)
    }else if(type == "loglog"){
        ## exp(-exp(log(-log(a)) +/- b)) = exp(-exp(log(-log(a)))exp(+/- b)) = exp(-(-log(a))exp(+/- b)) = exp(log(a)exp(+/- b)) = a ^ exp(+/- b)
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        out$lower <- estimate ^ exp(+ quantileSe)
        out$upper <- estimate ^ exp(- quantileSe)
    }else if(type == "cloglog"){
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
        out$lower <- 1 - (1-estimate) ^ exp(- quantileSe)
        out$upper <- 1 - (1-estimate) ^ exp(+ quantileSe)
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x) 
        ## fisher inverse  : g(x) = ( 1 - exp(-2x) ) / ( 1 + exp(-2x) )           
        out$lower <- tanh(atanh(estimate) - quantileSe)
        out$upper <- tanh(atanh(estimate) + quantileSe)
    }

    ## restrict to [min.value,max.value]
    if(!is.null(min.value)){        
        if(format == "matrix"){
            ## to keep matrix format even when object$survival contains only one line
            out$lower[] <- apply(out$lower,2,pmax,min.value)
        }else{
            out$lower <- pmax(out$lower,min.value)
        }
    }
    if(!is.null(max.value)){
        if(format == "matrix"){
            ## to keep matrix format even when object$survival contains only one line
            out$upper[] <- apply(out$upper,2,pmin,max.value)
        }else{
            out$upper <- pmin(out$upper,max.value)
        }
    }

    ## export
    return(out)
}

## * transformP
##' @title Compute the p.value after a transformation
##' @description Compute the p.value after a transformation.
##'
##' @param estimate [numeric vector] the estimate value before transformation.
##' @param se [numeric vector] the standard error after transformation.
##' @param null [numeric vector] the value of the estimate (before transformation) under the null hypothesis.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @export
transformP <- function(estimate, se, null, type){
    ## compute confidence intervals
    if(type == "none"){
        statistic <- ( estimate-null )/se
    }else if(type == "log"){
        statistic <- ( log(estimate) - log(null) )/se
    }else if(type == "loglog"){
        statistic <- ( log(-log(estimate)) - log(-log(null)) )/se
    }else if(type == "cloglog"){
        statistic <- ( log(-log(1-estimate)) - log(-log(1-null)) )/se
    }else if(type == "atanh"){
        statistic <- ( atanh(estimate) - atanh(null) )/se
    }

    return(2*(1-pnorm(abs(statistic))))
}

######################################################################
### transform.R ends here
