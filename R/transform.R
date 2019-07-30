### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 30 2018 (15:58) 
## Version: 
## Last-Updated: jul  2 2019 (11:05) 
##           By: Brice Ozenne
##     Update #: 176
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * transformSE
##' @title Compute Standard Errors after Transformation
##' @description Compute standard errors after transformation
##' based on the standard error before transformation.
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' 
##' @details Use a delta method to find the standard error after transformation.
##'
##' \code{se} and \code{estimate} must have same dimensions.
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
        newse <- se/estimate
    }else if(type == "loglog"){
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        newse <- se / ( - estimate * log(estimate) )
    }else if(type == "cloglog"){
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
        newse <- se / ( - (1-estimate) * log(1-estimate) )
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
        ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
        ##                         = dx/(1-x^2)
        ##               Var(f(x)) = Var(x)/(1-x^2)
        newse <- se / (1-estimate^2)
    }

    index0 <- which(se==0)
    if(length(index0)>0){
        newse[index0] <- 0
    }
    return(newse)
}

## * transformIID
##' @title Compute Influence Functions after Transformation
##' @description Compute influence functions after transformation
##' based on the influence function before transformation.
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param iid [numeric array] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' 
##' @details Use a delta method to find the standard error after transformation. \cr \cr
##'
##' The iid decomposition must contain have dimension [n.prediction,time,n.obs] and estimate [n.prediction,time].
##'
##' @export
transformIID <- function(estimate, iid, type){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.
    if(type == "none"){
        ## no change
        return(iid)
    }else if(type == "log"){
        ## formula 4.10 p 58 (Beyersmann et al. 2012)
        newiid <- sliceScale_cpp(iid, M = estimate)
    }else if(type == "loglog"){
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        newiid <- sliceScale_cpp(iid, M = - estimate * log(estimate) )
    }else if(type == "cloglog"){
        newiid <- sliceScale_cpp(iid, M = - (1 - estimate) * log(1-estimate) )
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
        ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
        ##                         = dx/(1-x^2)
        ##               Var(f(x)) = Var(x)/(1-x^2)
        newiid <- sliceScale_cpp(iid, M = 1 - estimate^2 )
    }
    index0 <- which(iid==0)
    if(length(index0)>0){
        newiid[index0] <- 0
    }
    return(newiid)
}

## * transformCI
##' @title Compute Confidence Intervals using a transformation
##' @description Compute confidence intervals using a transformation.
##' The resulting confidence interval is returned on the original case (i.e. back-transformed).
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error after transformation.
##' @param quantile [numeric vector] quantile that will be multiplied to each column of \code{se}.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @param min.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{min},
##' it will be set at \code{min}. 
##' @param max.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{max},
##' it will be set at \code{max}.
##'
##' @details \code{se} and \code{estimate} must have same dimensions.
##' 
##' @export
transformCI <- function(estimate, se, quantile, type, min.value, max.value){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.
    out <- list()
    quantileSe <- colMultiply_cpp(se, scale = quantile)
    
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
        ## to keep matrix format even when object$survival contains only one line
        out$lower[] <- apply(out$lower,2,pmax,min.value)
    }
    if(!is.null(max.value)){
        ## to keep matrix format even when object$survival contains only one line
        out$upper[] <- apply(out$upper,2,pmin,max.value)
    }

    ## export
    return(out)
}

## * transformP
##' @title Compute P-values After a Transformation
##' @description Compute the p-values after a transformation.
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error after transformation.
##' @param null [numeric] the value of the estimate (before transformation) under the null hypothesis.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' 
##' @details \code{se} and \code{estimate} must have same dimensions.
##' 
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

## * transformCIBP
##' @title Compute Confidence Intervals/Bands and P-values After a Transformation
##' @description Compute confidence intervals/bands and p-values after a transformation
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error before transformation.
##' @param iid [numeric array] the iid decomposition before transformation.
##' @param null [numeric] the value of the estimate (before transformation) under the null hypothesis.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @param conf.level [numeric, 0-1] Level of confidence.
##' @param nsim.band [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set before performing simulations for the confidence bands.
##' @param min.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{min},
##' it will be set at \code{min}. 
##' @param max.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{max},
##' it will be set at \code{max}.
##' @param ci [logical] should confidence intervals be computed.
##' @param band [logical] should confidence bands be computed.
##' @param p.value [logical] should p-values be computed.
##'
##' The iid decomposition must have dimensions [n.prediction,time,n.obs]
##' while estimate and se must have dimensions [n.prediction,time].
##' 
##' @export
transformCIBP <- function(estimate, se, iid, null,
                          conf.level, nsim.band, seed,
                          type, min.value, max.value, 
                          ci, band, p.value){

    out <- list()

    ## ** transformation
    if(type != "none"){
        ## standard error
        se <- transformSE(estimate = estimate,
                          se = se,
                          type = type)

        ## influence function
        if(band){
            iid <- transformIID(estimate = estimate,
                                iid = iid,
                                type = type)
        }
    }
    ## out$se <- se

    ## ** confidence intervals
    if(ci){
        zval <- stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1)

        out[c("lower","upper")] <- transformCI(estimate = estimate,
                                               se = se,
                                               quantile = rep(zval, NROW(estimate)),
                                               type = type,
                                               min.value = min.value,
                                               max.value = max.value)

    }

    ## ** confidence bands
    if(band[[1]] && nsim.band[[1]] > 0){

        ## find quantiles for the bands
        if(!is.na(seed)){set.seed(seed)}
        ## sqrt(rowSums(iid[1,,]^2))
        ## se[1,]
        index.firstSE <- which(colSums(se!=0)>0)[1]
        if(all(is.na(index.firstSE))){
            out$quantileBand <- rep(NA, NROW(se))
        }else if(index.firstSE!=1){
            ## sqrt(apply(iid[,index.firstSE:NCOL(se),,drop=FALSE]^2,1:2,sum))
            out$quantileBand <- confBandCox(iid = iid[,index.firstSE:NCOL(se),,drop=FALSE],
                                            se = se[,index.firstSE:NCOL(se),drop=FALSE],
                                            n.sim = nsim.band,
                                            conf.level = conf.level)
        }else{
            out$quantileBand <- confBandCox(iid = iid,
                                            se = se,
                                            n.sim = nsim.band,
                                            conf.level = conf.level)
        }
        out[c("lowerBand","upperBand")] <- transformCI(estimate = estimate,
                                                       se = se,
                                                       quantile = out$quantileBand,
                                                       type = type,
                                                       min.value = min.value,
                                                       max.value = max.value)

    }

    ## ** p.value
    if(p.value){
        out[["p.value"]] <- transformP(estimate = estimate,
                                       se = se,
                                       null = null,
                                       type = type)
    }
    
    ## ** check NA
    indexNA <- union(
        union(which(is.na(estimate)),
              which(is.nan(estimate))),
        union(which(is.na(se)),
              which(is.nan(se)))
    )
    
    if(length(indexNA)>0){
        if(ci){
            out$lower[indexNA] <- NA
            out$upper[indexNA] <- NA
        }
        if(band){ ## if cannot compute se at one time then remove confidence band at all times
            indexNA2 <- union(
                union(which(rowSums(is.na(estimate))>0),
                      which(rowSums(is.nan(estimate))>0)),
                union(which(rowSums(is.na(se))>0),
                      which(rowSums(is.nan(se))>0))
            )
            out$quantileBand[indexNA2] <- NA
            out$lowerBand[indexNA2,] <- NA
            out$upperBand[indexNA2,] <- NA
        }
    }

    ## ** export
    return(out)
}

######################################################################
### transform.R ends here
