### transform.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 30 2018 (15:58) 
## Version: 
## Last-Updated: sep 11 2024 (18:33) 
##           By: Brice Ozenne
##     Update #: 539
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
##' @noRd
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}, \code{"atanh"} (Fisher transform), or \code{"atanh2"} (modified Fisher transform for [0-1] variable).
##' 
##' @details Use a delta method to find the standard error after transformation.
##'
##' \code{se} and \code{estimate} must have same dimensions.
##'
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
    }else if(type == "atanh2"){
        newse <- 2*se / (1-4*(estimate-1/2)^2)
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
##' @noRd
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param iid [numeric array] the standard error before transformation.
##' @param type [character] the transforamtion.
##' Can be  \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}, \code{"atanh"} (Fisher transform), or \code{"atanh2"} (modified Fisher transform for [0-1] variable).
##' 
##' @details Use a delta method to find the standard error after transformation.
##'
##' The iid decomposition must contain have dimension [n.obs,time,n.prediction] and estimate [n.prediction,time].
##'
transformIID <- function(estimate, iid, type){
    ## Reference for the formula
    ## Beyersmann, Jan and Allignol, Arthur and Schumacher, Martin. Competing Risks and Multistate Models with R.
    ## Use R! Springer. 2012.
    # estimate has different time points in columns
    # and different covariate values in rows
    if(type == "none"){
        ## no change
        return(iid)
    }else{
        newiid <- aperm(iid, c(3,2,1))
        if(type == "log"){
            ## formula 4.10 p 58 (Beyersmann et al. 2012)
            ## newiid <- sliceScale_cpp(newiid, M = estimate)
            newiid <- sweep(newiid,MARGIN = 1:2,FUN = "/",STATS = estimate)
        }else if(type == "loglog"){
            ## formula 4.16 p 59 (Beyersmann et al. 2012)
            ## newiid <- sliceScale_cpp(newiid, M = - estimate * log(estimate) )
            newiid <- sweep(newiid,MARGIN = 1:2,FUN = "/",STATS = - estimate * log(estimate))
        }else if(type == "cloglog"){
            ## newiid <- sliceScale_cpp(newiid, M = - (1 - estimate) * log(1-estimate) )
            newiid <- sweep(newiid,MARGIN = 1:2,FUN = "/",STATS = - (1 - estimate) * log(1-estimate))
            ## formula 4.21 p 62 (Beyersmann et al. 2012)
        }else if(type == "atanh"){
            ## fisher transform: f(x) = 1/2 log(1+x) - 1/2 log(1-x)
            ##                   df(x) = dx/(2+2x) + dx/(2-2x) = dx(1/(1+x)+1/(1-x))/2
            ##                         = dx/(1-x^2)
            ##               Var(f(x)) = Var(x)/(1-x^2)
            ## newiid <- sliceScale_cpp(newiid, M = 1 - estimate^2 )
            newiid <- sweep(newiid,MARGIN = 1:2,FUN = "/",STATS = 1 - estimate^2)
        }else if(type == "atanh2"){
            ## newiid <- sliceScale_cpp(2*newiid, M = 1 - 4*(estimate-1/2)^2 )
            newiid <- sweep(2*newiid,MARGIN = 1:2,FUN = "/",STATS = 1 - 4*(estimate-1/2)^2 )
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
##' @noRd
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error after transformation.
##' @param quantile [numeric vector] quantile that will be multiplied to each column of \code{se}.
##' @param type [character] the transforamtion.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}, \code{"atanh"} (Fisher transform), or \code{"atanh2"} (modified Fisher transform for [0-1] variable).
##' @param min.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{min},
##' it will be set at \code{min}. 
##' @param max.value [numeric] if not \code{NULL} and the lower bound of the confidence interval is below \code{max},
##' it will be set at \code{max}.
##'
##' @details \code{se} and \code{estimate} must have same dimensions.
##' 
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
    }else if(type == "atanh2"){
        out$lower <- tanh(atanh(2*(estimate-0.5)) + quantileLowerSe)/2 + 1/2
        out$upper <- tanh(atanh(2*(estimate-0.5)) + quantileUpperSe)/2 + 1/2
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
##' @noRd
##'
##' @param estimate [numeric matrix] the estimate value before transformation.
##' @param se [numeric matrix] the standard error after transformation.
##' @param null [numeric] the value of the estimate (before transformation) under the null hypothesis.
##' @param type [character] the transforamtion.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}, \code{"atanh"} (Fisher transform), or \code{"atanh2"} (modified Fisher transform for [0-1] variable).
##' @param alternative [character] a character string specifying the alternative hypothesis,
##' must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
##' 
##' @details \code{se} and \code{estimate} must have same dimensions.
##' 
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
    }else if(type == "atanh2"){
        statistic <- ( atanh(2*(estimate - null)) )/se
    }
    statistic[estimate==null] <- 0 ## deal with 0 estimate with 0 variance

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
##' Can be \code{"log"}, \code{"loglog"}, \code{"cloglog"}, or \code{"atanh"} (Fisher transform), or \code{"atanh2"} (modified Fisher transform for [0-1] variable).
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
##' @param df [integer, >0] optional. Degrees of freedom used for the student distribution of the test statistic. If not specified, use a normal distribution instead.
##' 
##' @details The iid decomposition must have dimensions [n.obs,time,n.prediction]
##' while estimate and se must have dimensions [n.prediction,time].
##'
##' Single step max adjustment for multiple comparisons, i.e. accounting for the correlation between the test statistics but not for the ordering of the tests, can be performed setting the arguemnt \code{method.band} to \code{"maxT-integration"} or \code{"maxT-simulation"}. The former uses numerical integration (\code{pmvnorm} and \code{qmvnorm} to perform the adjustment while the latter using simulation. Both assume that the test statistics are jointly normally distributed. 
##'
##' @examples
##' set.seed(10)
##' n <- 100
##' X <- rnorm(n)
##' 
##' res2sided <- transformCIBP(estimate = mean(X), se = cbind(sd(X)/sqrt(n)), null = 0,
##'               type = "none", ci = TRUE, conf.level = 0.95, alternative = "two.sided",
##'               min.value = NULL, max.value = NULL, band = FALSE,
##'               p.value = TRUE, seed = 10, df = n-1)
##' 
##' resLess <- transformCIBP(estimate = mean(X), se = cbind(sd(X)/sqrt(n)), null = 0,
##'               type = "none", ci = TRUE, conf.level = 0.95, alternative = "less",
##'               min.value = NULL, max.value = NULL, band = FALSE,
##'               p.value = TRUE, seed = 10, df = n-1)
##' 
##' resGreater <- transformCIBP(estimate = mean(X), se = cbind(sd(X)/sqrt(n)), null = 0,
##'               type = "none", ci = TRUE, conf.level = 0.95, alternative = "greater",
##'               min.value = NULL, max.value = NULL, band = FALSE,
##'               p.value = TRUE, seed = 10, df = n-1)
##'
##'
##' ## comparison with t-test
##' GS <- t.test(X, alternative = "two.sided")
##' res2sided$p.value - GS$p.value
##' unlist(res2sided[c("lower","upper")]) - GS$conf.int
##' 
##' GS <- t.test(X, alternative = "less")
##' resLess$p.value - GS$p.value
##' unlist(resLess[c("lower","upper")]) - GS$conf.int
##' 
##' GS <- t.test(X, alternative = "greater")
##' resGreater$p.value - GS$p.value
##' unlist(resGreater[c("lower","upper")]) - GS$conf.int
##'
##' @export
transformCIBP <- function(estimate, se, iid, null,
                          conf.level, alternative,                        
                          ci, type, min.value, max.value,
                          band, method.band, n.sim, seed,
                          p.value, df = NULL){

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
    if(band){
        method.band <- match.arg(method.band, choices = c(setdiff(p.adjust.methods,"none"),"maxT-integration","maxT-simulation"))
        if(all((abs(se[!is.na(se)])<1e-12))){
            method.band <- "bonferroni"
        }
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
    if(band>0){
        if(method.band %in% c("maxT-integration","maxT-simulation")){
            iid <- transformIID(estimate = estimate,
                                iid = iid,
                                type = type)
        }
    }
    
    ## ** normalize influence function and statistic
    if(band){
        if(method.band %in% c("maxT-integration","maxT-simulation")){
            n.sample <- dim(iid)[1]
            n.time <- dim(iid)[2]

            ## times with 0 variance (to be removed in further calculation as they introduce singularities)
            index.keep <- which(colSums(abs(se)>1e-12, na.rm = TRUE)>0)[1]:utils::tail(which(colSums(abs(se)>1e-12, na.rm = TRUE)>0),1)
            iid.norm <- array(NA, dim = c(length(index.keep), dim(iid)[1], dim(iid)[3]))
            for(iC in 1:n.contrast){ ## iC <- 1
                if(length(index.keep)==1){
                    iid.norm[,,iC] <- t(iid[,index.keep,iC] / se[iC,index.keep])
                }else{
                    iid.norm[,,iC] <- t(rowScale_cpp(iid[,index.keep,iC], scale = se[iC,index.keep]))
                }
            }
            if(band==1){
                if(n.time==1){
                    rho <- lapply(diag(crossprod(iid.norm[1,,])),as.matrix)
                }else{
                    rho <- lapply(1:n.contrast, function(iC){tcrossprod(iid.norm[,,iC])})
                }
            }else if(band == 2){
                if(n.time==1){
                    rho <- crossprod(iid.norm[1,,])
                }else{
                    rho <- tcrossprod(do.call(rbind,lapply(1:n.contrast, function(iC){iid.norm[,,iC]})))
                }
            }
        }
    }

    ## ** unadjusted: confidence intervals and p-value
    if(ci){
        if(alternative == "two.sided"){
            if(is.null(df)){
                zval <- c(stats::qnorm((1-conf.level)/2, mean = 0, sd = 1),
                          stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1))
            }else{
                zval <- c(stats::qt((1-conf.level)/2, df = df),
                          stats::qt(1 - (1-conf.level)/2, df = df))
            }
        }else if(alternative == "greater"){
            if(is.null(df)){
                zval <- c(stats::qnorm(1-conf.level, mean = 0, sd = 1),
                          Inf)
            }else{
                zval <- c(stats::qt(1-conf.level, df = df),
                          Inf)
            }
        }else if(alternative == "less"){
            if(is.null(df)){
                zval <- c(-Inf,
                          stats::qnorm(conf.level, mean = 0, sd = 1))
            }else{
                zval <- c(-Inf,
                          stats::qt(conf.level, df = df))
            }
                      
        }
        out[c("lower","upper")] <- transformCI(estimate = estimate,
                                               se = se,
                                               quantile = matrix(zval, nrow = n.contrast, ncol = 2, byrow = TRUE),
                                               type = type,
                                               min.value = min.value,
                                               max.value = max.value)

        if(p.value){
            if(alternative == "two.sided"){
                if(is.null(df)){
                    out[["p.value"]] <- 2*(1-stats::pnorm(abs(statistic)))
                }else{
                    out[["p.value"]] <- 2*(1-stats::pt(abs(statistic), df = df))
                }
            }else if(alternative == "less"){
                if(is.null(df)){
                    out[["p.value"]] <- stats::pnorm(statistic)
                }else{
                    out[["p.value"]] <- stats::pt(statistic, df = df)
                }
            }else if(alternative == "greater"){
                if(is.null(df)){
                    out[["p.value"]] <- 1-stats::pnorm(statistic)
                }else{
                    out[["p.value"]] <- 1-stats::pt(statistic, df = df)
                }
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

                        resQ <- mvtnorm::qmvnorm(p = conf.level, mean = rep(0, NCOL(rho[[iC]])),
                                                 sigma = rho[[iC]], tail = switch(alternative,
                                                                                  "two.sided" = "both.tails",
                                                                                  "less" = "upper.tail", ## 'upper.tail' gives x with P[X > x] = p = P[x < X < Inf]
                                                                                  "greater" = "lower.tail"))$quantile
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
                    if(any(is.na(rho)) && any(se==0)){
                        rho <- rho[se>0,se>0,drop=FALSE]
                    }
                    resQ <- mvtnorm::qmvnorm(p = conf.level, mean = rep(0, NCOL(rho)),
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
                                              iid = iid.norm,
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
            out$adj.p.value <- matrix(NA, nrow = n.contrast, ncol = NCOL(estimate))

            if(method.band %in% p.adjust.methods){
                if(band == 1){
                    for(iC in 1:NROW(out[["p.value"]])){
                        out[["adj.p.value"]][iC,] <- stats::p.adjust(out[["p.value"]][iC,], method = method.band)
                    }
                }else if(band == 2){
                    out[["adj.p.value"]] <- stats::p.adjust(out[["p.value"]], method = method.band)
                }
            }else if(method.band == "maxT-integration"){

                if(band == 1){
                    for(iC in 1:n.contrast){ ## iC <- 2
                        iN.test <- NCOL(rho[[iC]])

                        out$adj.p.value[iC,index.keep] <-  sapply(statistic[iC,index.keep], function(iT){ ## iT
                            if(alternative=="two.sided"){
                                return(1-mvtnorm::pmvnorm(lower=rep(-abs(iT), iN.test), upper=rep(abs(iT), iN.test),
                                                          mean=rep(0, iN.test), sigma = rho[[iC]]))
                            }else if(alternative=="greater"){
                                return(1-mvtnorm::pmvnorm(lower=rep(-Inf, iN.test), upper=rep(iT, iN.test),
                                                          mean=rep(0, iN.test), sigma = rho[[iC]]))
                            }else if(alternative=="less"){
                                return(1-mvtnorm::pmvnorm(lower=rep(iT, iN.test), upper=rep(Inf, iN.test),
                                                          mean=rep(0, iN.test), sigma = rho[[iC]]))
                            }
                        })

                        if(length(index.keep)>0 && all(stats::na.omit(out$p.value[iC,-index.keep])==1)){
                            out$adj.p.value[iC,-index.keep] <- out$p.value[iC,-index.keep]
                        }
                    }
                }else if(band == 2){
                    for(iC in 1:n.contrast){ ## iC <- 2
                        iN.test <- NCOL(rho[[iC]])

                        out$adj.p.value[iC,index.keep] <- sapply(statistic[iC,index.keep], function(iT){
                            if(alternative=="two.sided"){
                                return(1-mvtnorm::pmvnorm(lower=rep(-abs(iT), iN.test), upper=rep(abs(iT), iN.test),
                                                          mean=rep(0, iN.test), sigma = rho))
                            }else if(alternative=="greater"){
                                return(1-mvtnorm::pmvnorm(lower=rep(-Inf, iN.test), upper=rep(iT,n.test),
                                                          mean=rep(0, iN.test), sigma = rho))
                            }else if(alternative=="less"){
                                return(1-mvtnorm::pmvnorm(lower=rep(iT, iN.test), upper=rep(Inf, iN.test),
                                                          mean=rep(0, iN.test), sigma = rho))
                            }
                        })

                        if(length(index.keep)>0 && all(stats::na.omit(out$p.value[iC,-index.keep])==1)){
                            out$adj.p.value[iC,-index.keep] <- out$p.value[iC,-index.keep]
                        }
                    }
                }
            }else if(method.band == "maxT-simulation"){
                out$adj.p.value[,index.keep] <- pProcess_cpp(nSample = n.sample,
                                                             nContrast = n.contrast,
                                                             nTime = length(index.keep),
                                                             nSim = n.sim,
                                                             value = statistic[,index.keep,drop=FALSE],
                                                             iid =  iid.norm,
                                                             alternative = switch(alternative,
                                                                                  "two.sided" = 3,
                                                                                  "greater" = 2,
                                                                                  "less" = 1),
                                                             global = (band == 2)                             
                                                             )

                if(length(index.keep)>0 && all(stats::na.omit(out$p.value[,-index.keep])==1)){
                    out$adj.p.value[,-index.keep] <- out$p.value[,-index.keep]
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
        if(band){
            out$lowerBand[indexNA] <- NA
            out$upperBand[indexNA] <- NA
            if(p.value){
                out$adj.p.value[indexNA] <- NA
            }
            ## if cannot compute se at one time then remove confidence band at all times
            ## indexNA2 <- union(
            ##     union(which(rowSums(is.na(estimate))>0),
            ##           which(rowSums(is.nan(estimate))>0)),
            ##     union(which(rowSums(is.na(se))>0),
            ##           which(rowSums(is.nan(se))>0))
            ## )
            ## out$quantileBand[indexNA2] <- NA
            ## out$lowerBand[indexNA2,] <- NA
            ## out$upperBand[indexNA2,] <- NA
            ## if(p.value){
            ##     out$adj.p.value[indexNA2,] <- NA
            ## }
        }
    }

    ## ** export
    class(out) <- "transformCIBP"
    return(out)
}

## * as.data.table.transformCIBP
##' @export
as.data.table.transformCIBP <- function(x, keep.rownames = FALSE, ...){

    ## add extra argument (should be of the correct size)
    dots <- list(...)
    x <- c(x,dots)

    ## prepare
    name.x <- names(x)
    nameRow.x <- setdiff(names(x),"quantileBand")

    n.endpoint <- NROW(x[[nameRow.x[1]]])
    n.times <- NCOL(x[[nameRow.x[1]]])
    if(any(sapply(dots,is.matrix)==FALSE)){
        stop("Extra arguments must be matrices \n")
    }
    if(any(sapply(dots,NROW)!=n.endpoint)){
        stop("Extra arguments must have ",n.endpoint," row(s) \n")
    }
    if(any(sapply(dots,NCOL)!=n.times)){
        stop("Extra arguments must have ",n.times," column(s) \n")
    }
    
    ## merge
    dt <- NULL
    for(iEndpoint in 1:n.endpoint){ ## iEndpoint <- 1
        iDT <- as.data.table(lapply(x[nameRow.x], function(iE){iE[iEndpoint,]}))
        if("quantileBand" %in% name.x){
            iDT[,c("quantileBand") := x$quantileBand[iEndpoint]]
        }
        iDT[,c("row") := iEndpoint]
        iDT[,c("time") := 1:.N]
        dt <- rbind(dt, iDT)
    }
    setcolorder(dt, neworder = c("row","time",name.x))
    return(dt)
}

######################################################################
### transform.R ends here
