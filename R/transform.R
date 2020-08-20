### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 30 2018 (15:58) 
## Version: 
## Last-Updated: aug 20 2020 (15:49) 
##           By: Brice Ozenne
##     Update #: 394
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
##' @details Use a delta method to find the standard error after transformation.
##'
##' The iid decomposition must contain have dimension [n.obs,time,n.prediction] and estimate [n.prediction,time].
##'
##' @export
transformIID <- function(estimate, iid, type){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.
    if(type == "none"){
        ## no change
        return(iid)
    }else{
        newiid <- aperm(iid, c(3,2,1))
        if(type == "log"){
            ## formula 4.10 p 58 (Beyersmann et al. 2012)
            newiid <- sliceScale_cpp(newiid, M = estimate)
        }else if(type == "loglog"){
            ## formula 4.16 p 59 (Beyersmann et al. 2012)
            newiid <- sliceScale_cpp(newiid, M = - estimate * log(estimate) )
        }else if(type == "cloglog"){
            newiid <- sliceScale_cpp(newiid, M = - (1 - estimate) * log(1-estimate) )
            ## formula 4.21 p 62 (Beyersmann et al. 2012)
        }else if(type == "atanh"){
            ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
            ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
            ##                         = dx/(1-x^2)
            ##               Var(f(x)) = Var(x)/(1-x^2)
            newiid <- sliceScale_cpp(newiid, M = 1 - estimate^2 )
        }
        newiid <- aperm(newiid, c(3,2,1))
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
    quantileLowerSe <- colMultiply_cpp(se, scale = quantile[,1])
    quantileUpperSe <- colMultiply_cpp(se, scale = quantile[,2])
    
    ## compute confidence intervals
    if(type == "none"){
        out$lower <- estimate + quantileLowerSe
        out$upper <- estimate + quantileUpperSe
    }else if(type == "log"){
        ## a * exp +/-b = exp(log(a) +/- b)
        ## formula 4.10 p 58 (Beyersmann et al. 2012)
        out$lower <- estimate * exp(quantileLowerSe)
        out$upper <- estimate * exp(quantileUpperSe)
    }else if(type == "loglog"){
        ## exp(-exp(log(-log(a)) +/- b)) = exp(-exp(log(-log(a)))exp(+/- b)) = exp(-(-log(a))exp(+/- b)) = exp(log(a)exp(+/- b)) = a ^ exp(+/- b)
        ## formula 4.16 p 59 (Beyersmann et al. 2012)
        out$lower <- estimate ^ exp(- quantileLowerSe)
        out$upper <- estimate ^ exp(- quantileUpperSe)
    }else if(type == "cloglog"){
        ## formula 4.21 p 62 (Beyersmann et al. 2012)
        out$lower <- 1 - (1-estimate) ^ exp(quantileLowerSe)
        out$upper <- 1 - (1-estimate) ^ exp(quantileUpperSe)
    }else if(type == "atanh"){
        ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x) 
        ## fisher inverse  : g(x) = ( 1 - exp(-2x) ) / ( 1 + exp(-2x) )           
        out$lower <- tanh(atanh(estimate) + quantileLowerSe)
        out$upper <- tanh(atanh(estimate) + quantileUpperSe)
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

## * transformT
##' @title Compute T-statistic after a Transformation
##' @description Compute T-statistic after a transformation.
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error after transformation.
##' @param null [numeric] the value of the estimate (before transformation) under the null hypothesis.
##' @param type [character] the transforamtion.
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform).
##' @param alternative [character] a character string specifying the alternative hypothesis,
##' must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
##' 
##' @details \code{se} and \code{estimate} must have same dimensions.
##' 
##' @export
transformT <- function(estimate, se, null, type, alternative){
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

    return(statistic)
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
##' @param alternative [character] a character string specifying the alternative hypothesis,
##' must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
##' @param ci [logical] should confidence intervals be computed.
##' @param min.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{min},
##' it will be set at \code{min}. 
##' @param max.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{max},
##' it will be set at \code{max}.
##' @param band [integer 0,1,2] When non-0, the confidence bands are computed for each contrasts (\code{band=1}) or over all contrasts (\code{band=2}). 
##' @param method.band [character] method used to adjust for multiple comparisons.
##' Can be any element of \code{p.adjust.methods} (e.g. \code{"holm"}), \code{"maxT-integration"}, or \code{"maxT-simulation"}. 
##' @param n.sim [integer, >0] the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed [integer, >0] seed number set before performing simulations for the confidence bands.
##' @param p.value [logical] should p-values and adjusted p-values be computed. Only active if \code{ci=TRUE} or \code{band>0}.
##' 
##' @details The iid decomposition must have dimensions [n.obs,time,n.prediction]
##' while estimate and se must have dimensions [n.prediction,time].
##'
##' Single step max adjustment for multiple comparisons, i.e. accounting for the correlation between the test statistics but not for the ordering of the tests, can be performed setting the arguemnt \code{method.band} to \code{"maxT-integration"} or \code{"maxT-simulation"}. The former uses numerical integration (\code{pmvnorm} and \code{qmvnorm} to perform the adjustment while the latter using simulation. Both assume that the test statistics are jointly normally distributed. 
##' 
##' @export
transformCIBP <- function(estimate, se, iid, null,
                          conf.level, alternative,                        
                          ci, type, min.value, max.value,
                          band, method.band, n.sim, seed,
                          p.value){

    p.adjust.methods <-  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    out <- list()
    if(band %in% c(0,1,2) == FALSE){
        stop("Incorrect value for argument \'band\'. Should be 0, 1, or 2. \n")
    }
    if(band == 1){
        n.test <- NCOL(estimate)
    }else if(band == 2){
        n.test <- length(estimate)
    }
    n.contrast <- NROW(estimate)
    method.band <- match.arg(method.band, choices = c(setdiff(p.adjust.methods,"none"),"maxT-integration","maxT-simulation"))
    if(all((abs(se[!is.na(se)])<1e-12))){
        method.band <- "bonferroni"
    }
    if(!is.na(seed)){set.seed(seed)}
    alternative <- match.arg(alternative, choices = c("two.sided","less","greater"))
    
    ## ** transformation
    ## standard error
    se <- transformSE(estimate = estimate,
                      se = se,
                      type = type)

    ## test statistic
    if(p.value){
        statistic <- transformT(estimate = estimate,
                                se = se,
                                null = null,
                                type = type)
    }
    
    ## influence function
    if(band>0 && method.band %in% c("maxT-integration","maxT-simulation")){
        iid <- transformIID(estimate = estimate,
                            iid = iid,
                            type = type)
    }
    
    ## ** normalize influence function and statistic
    if(band>0 && method.band %in% c("maxT-integration","maxT-simulation")){
        n.sample <- dim(iid)[1]
        n.time <- dim(iid)[2]

        ## times with 0 variance (to be removed in further calculaltion as they introduce singularities)
        index.keep <- which(colSums(abs(se)>1e-12, na.rm = TRUE)>0)[1]:NCOL(se)
        iid.norm <- array(NA, dim = c(dim(iid)[1],length(index.keep), dim(iid)[3]))
        for(iC in 1:n.contrast){ ## iC <- 1
            if(length(index.keep)==1){
                iid.norm[,,iC] <- iid[,index.keep,iC] / se[iC,index.keep]
            }else{
                iid.norm[,,iC] <- rowScale_cpp(iid[,index.keep,iC], scale = se[iC,index.keep])
            }
        }
        if(band==1){
            rho <- lapply(1:n.contrast, function(iC){crossprod(iid.norm[,,iC])})
        }else if(band == 2){
            rho <- crossprod(do.call(cbind,lapply(1:n.contrast, function(iC){iid.norm[,,iC]})))
        }
    }

    ## ** unadjusted: confidence intervals and p-value
    if(ci){
        if(alternative == "two.sided"){
            zval <- c(stats::qnorm((1-conf.level)/2, mean = 0, sd = 1),
                      stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1))
        }else if(alternative == "greater"){
            zval <- c(-Inf,
                      stats::qnorm(conf.level, mean = 0, sd = 1))
                      
        }else if(alternative == "less"){
            zval <- c(stats::qnorm(1-conf.level, mean = 0, sd = 1),
                      Inf)
                      
        }
        out[c("lower","upper")] <- transformCI(estimate = estimate,
                                               se = se,
                                               quantile = matrix(zval, nrow = n.contrast, ncol = 2, byrow = TRUE),
                                               type = type,
                                               min.value = min.value,
                                               max.value = max.value)
        if(p.value){
            if(alternative == "two.sided"){
                out[["p.value"]] <- 2*(1-pnorm(abs(statistic)))
            }else if(alternative == "less"){
                out[["p.value"]] <-- pnorm(statistic)
            }else if(alternative == "greater"){
                out[["p.value"]] <- -pnorm(statistic)
            }
        }
    }

    ## ** adjusted: confidence bands and adjusted p.value
    if(band>0){
        
        if(method.band == "bonferroni"){
            if(alternative == "two.sided"){
                quantileBand <- c(stats::qnorm((1-conf.level)/(2*n.test), mean = 0, sd = 1),
                                  stats::qnorm(1 - (1-conf.level)/(2*n.test), mean = 0, sd = 1))
            }else if(alternative == "greater"){
                quantileBand <- c(-Inf,
                                  stats::qnorm(1-(1-conf.level)/n.test, mean = 0, sd = 1))
                      
            }else if(alternative == "less"){
                quantileBand <- c(stats::qnorm((1-conf.level)/n.test, mean = 0, sd = 1),
                                  Inf)
                      
            }
            quantileBand <- matrix(quantileBand, byrow = TRUE, nrow = n.contrast, ncol = 2)
        }else if(method.band %in% c("maxT-integration","maxT-simulation")){
            quantileBand <- matrix(NA, nrow = n.contrast, ncol = 2)
            if(method.band == "maxT-integration"){
                if(band == 1){
                    for(iC in 1:n.contrast){ ## iC <- 1
                        if(n.test==1){
                            resQ <- mvtnorm::qmvnorm(p = conf.level, mean = rep(0, n.test),
                                                     sigma = rho[[iC]], tail = switch(alternative,
                                                                                      "two.sided" = "both.tails",
                                                                                      "less" = "upper.tail", ## 'upper.tail' gives x with P[X > x] = p = P[x < X < Inf]
                                                                                      "greater" = "lower.tail"))$quantile
                        }else{
                            resQ <- mvtnorm::qmvnorm(p = conf.level, mean = rep(0, n.test),
                                                     cor = rho[[iC]], tail = switch(alternative,
                                                                                    "two.sided" = "both.tails",
                                                                                    "less" = "upper.tail", ## 'upper.tail' gives x with P[X > x] = p = P[x < X < Inf]
                                                                                    "greater" = "lower.tail"))$quantile
                        }
                        if(alternative == "two.sided"){
                            quantileBand[iC,1] <- -resQ
                            quantileBand[iC,2] <- resQ
                        }else if(alternative == "greater"){
                            quantileBand[iC,1] <- -Inf
                            quantileBand[iC,2] <- resQ
                        }else if(alternative == "less"){
                            quantileBand[iC,1] <- resQ
                            quantileBand[iC,2] <- Inf
                        }
                    }
                }else if(band == 2){
                    resQ <- mvtnorm::qmvnorm(p = conf.level, mean = rep(0, n.test),
                                             cor = rho, tail = switch(alternative,
                                                                      "two.sided" = "both.tails",
                                                                      "less" = "upper.tail", ## 'upper.tail' gives x with P[X > x] = p = P[x < X < Inf]
                                                                      "greater" = "lower.tail"))$quantile
                    if(alternative == "two.sided"){
                        quantileBand[,1] <- -resQ
                        quantileBand[,2] <- resQ
                    }else if(alternative == "greater"){
                        quantileBand[,1] <- -Inf
                        quantileBand[,2] <- resQ
                    }else if(alternative == "less"){
                        quantileBand[,1] <- resQ
                        quantileBand[,2] <- Inf
                    }
                }
            
            }else if(method.band == "maxT-simulation"){
                resCpp <- quantileProcess_cpp(nSample = n.sample,
                                              nContrast = n.contrast,
                                              nSim = n.sim,
                                              iid = aperm(iid.norm, perm = c(2, 1, 3)),
                                              alternative = switch(alternative,
                                                                   "two.sided" = 3,
                                                                   "greater" = 2,
                                                                   "less" = 1),
                                              global = (band == 2),
                                              confLevel = conf.level)
                if(alternative == "two.sided"){
                    quantileBand[,1] <- -resCpp
                    quantileBand[,2] <- resCpp
                }else if(alternative == "greater"){
                    quantileBand[,1] <- -Inf
                    quantileBand[,2] <- resCpp
                }else if(alternative == "less"){
                    quantileBand[,1] <- resCpp
                    quantileBand[,2] <- Inf
                }
            }
            
        }

        if(method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            if(alternative == "less"){
                out$quantileBand <- -quantileBand[,1]
            }else{
                out$quantileBand <- quantileBand[,2]
            }
            out[c("lowerBand","upperBand")] <- transformCI(estimate = estimate,
                                                           se = se,
                                                           quantile = quantileBand,
                                                           type = type,
                                                           min.value = min.value,
                                                           max.value = max.value)
        }

        if(p.value){
            if(method.band %in% p.adjust.methods){
                if(band == 1){
                    out$adj.p.value <- matrix(NA, nrow = n.contrast, ncol = NCOL(estimate))
                    for(iC in 1:NROW(out[["p.value"]])){
                        out[["adj.p.value"]][iC,] <- stats::p.adjust(out[["p.value"]][iC,], method = method.band)
                    }
                }else if(band == 2){
                    out[["adj.p.value"]] <- stats::p.adjust(out[["p.value"]], method = method.band)
                }
            }else if(method.band == "maxT-integration"){
                out$adj.p.value <- matrix(NA, nrow = n.contrast, ncol = NCOL(estimate))

                if(band == 1){
                    for(iC in 1:n.contrast){ ## iC <- 1
                        out$adj.p.value[iC,] <- sapply(statistic[iC,], function(iT){
                            if(alternative=="two.sided"){
                                return(1-mvtnorm::pmvnorm(lower=-abs(iT), upper=abs(iT),
                                                 mean=rep(0, n.test), corr = rho[[iC]]))
                            }else if(alternative=="greater"){
                                return(mvtnorm::pmvnorm(lower=iT, upper=Inf,
                                               mean=rep(0, n.test), corr = rho[[iC]]))
                            }else if(alternative=="less"){
                                return(mvtnorm::pmvnorm(lower=-Inf, upper=iT,
                                               mean=rep(0, n.test), corr = rho[[iC]]))
                            }
                        })
                    }
                }else if(band == 2){
                    out$adj.p.value[iC,] <- sapply(statistic, function(iT){
                        if(alternative=="two.sided"){
                            return(1-mvtnorm::pmvnorm(lower=-abs(iT), upper=abs(iT),
                                             mean=rep(0, n.test), corr = rho))
                        }else if(alternative=="greater"){
                            return(mvtnorm::pmvnorm(lower=iT, upper=Inf,
                                           mean=rep(0, n.test), corr = rho))
                        }else if(alternative=="less"){
                            return(mvtnorm::pmvnorm(lower=-Inf, upper=iT,
                                           mean=rep(0, n.test), corr = rho))
                        }
                    })
                }
            }else if(method.band == "maxT-simulation"){
                out$adj.p.value <- matrix(NA, nrow = n.contrast, ncol = NCOL(estimate))
                for(iC in 1:n.contrast){ ## iC <- 2
                    if(band == 1){
                        iid.norm.tempo <- aperm(iid.norm[,,iC,drop=FALSE], c(2,1,3))
                    }else{
                        iid.norm.tempo <- aperm(iid.norm, c(2,1,3))
                    }
                    out$adj.p.value[iC,] <- sapply(statistic[iC,], function(iT){
                        mean(sampleMaxProcess_cpp(nSample = n.sample,
                                                  nContrast = dim(iid.norm.tempo)[3],
                                                  nSim = n.sim,
                                                  value = matrix(iT, nrow = n.time, ncol = dim(iid.norm.tempo)[3]),
                                                  iid =  iid.norm.tempo,
                                                  global = (band == 2),
                                                  type = 1,
                                                  alternative = switch(alternative,
                                                                       "two.sided" = 3,
                                                                       "greater" = 2,
                                                                       "less" = 1)
                                                  )>=0)
                    })
                }
            }
        }

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
            if(p.value){
                out$p.value[indexNA] <- NA
            }
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
            if(p.value){
                out$adj.p.value[indexNA2,] <- NA
            }
        }
    }

    ## ** export
    return(out)
}

######################################################################
### transform.R ends here
