### confint.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 24 2018 (17:53) 
##           By: Brice Ozenne
##     Update #: 145
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Confidence Intervals and Confidence Bands for the predicted Survival/Cumulative Hazard
##' @description Confidence intervals and confidence Bands for the predicted survival/cumulative Hazard.
##' @name confint.predictCox
##' 
##' @param object A \code{predictCox} object, i.e. output of the \code{predictCox} function.
##' @param conf.level Level of confidence.
##' @param type the type of predicted value for which the confidence intervals should be output.
##' Can be \code{"survival"} or \code{"cumhazard"}.
##' @param transform.cumhazard the transformation used to improve coverage of the confidence intervals for the cumlative hazard in small samples.
##' @param transform.survival the transformation used to improve coverage of the confidence intervals for the survival in small samples.
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
##' @export
confint.predictCox <- function(object,
                               conf.level = 0.95,
                               type = NULL,
                               nsim.band = 1e4,
                               transform.cumhazard = "log",
                               transform.survival = "cloglog",
                               seed = NA,
                               ...){

    ## check arguments
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
        if(!is.null(object$transform.cumhazard) && object$transform.cumhazard != "none"){
            stop("Cannot work with standard errors that have already been transformed \n")
        }
        object$transform.cumhazard <- match.arg(transform.cumhazard, c("none","log"))
    }
    if("survival" %in% type){
        if(!is.null(object$transform.survival) && object$transform.survival != "none"){
            stop("Cannot work with standard errors that have already been transformed \n")
        }
        object$transform.survival <- match.arg(transform.survival, c("none","log","loglog","cloglog"))
    }
    ## extract information
    se <- object$se
    band <- object$band

    ## quantile
    zval <- stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1)

    ## standard error and influence function
    if("cumhazard" %in% type){
        if(object$transform.cumhazard == "none"){
            ## no change
        }else if(object$transform.cumhazard == "log"){
            ## formula 4.10 p 58 (Beyersmann et al. 2012)
            object$cumhazard.se <- object$cumhazard.se/object$cumhazard
            if(band){
                object$cumhazard.iid <- sweep(object$cumhazard.iid, MARGIN = 1:2, FUN = "/", STATS = object$cumhazard)
            }
        }
    }
        
    if("survival" %in% type){
        if(object$transform.survival == "none"){
            ## no change
        }else if(object$transform.survival == "log"){
            ## formula 4.10 p 58 (Beyersmann et al. 2012)
            object$survival.se <- object$survival.se / object$survival
            if(band){
                object$survival.iid <- sweep(object$survival.iid, MARGIN = 1:2, FUN = "/", STATS = object$survival)
            }
        }else if(object$transform.survival == "loglog"){
            ## formula 4.16 p 59 (Beyersmann et al. 2012)
            object$survival.se <- object$survival.se / (object$survival * (-log(object$survival)))
            if(band){
                object$survival.iid <- sweep(object$survival.iid, MARGIN = 1:2, FUN = "/", STATS = object$survival * log(object$survival))
            }
        }else if(object$transform.survival == "cloglog"){
            ## formula 4.21 p 62 (Beyersmann et al. 2012)
            object$survival.se <- object$survival.se / ((1-object$survival) * (-log(1-object$survival)))
            if(band){
                object$survival.iid <- sweep(object$survival.iid, MARGIN = 1:2, FUN = "/", STATS = (1-object$survival) * log(1-object$survival))
            }
        }
    }

    
    ## confidence intervals
    if(object$se){
        if("cumhazard" %in% type){
            if(object$transform.cumhazard == "none"){
                object$cumhazard.lower <- matrix(NA, nrow = NROW(object$cumhazard), ncol = NCOL(object$cumhazard)) # to keep matrix format even when object$cumhazard contains only one line
                object$cumhazard.lower[] <- apply(object$cumhazard - zval * object$cumhazard.se,2,pmax,0) # to keep matrix format even when object$cumhazard contains only one line
                object$cumhazard.upper <- object$cumhazard + zval * object$cumhazard.se
            }else if(object$transform.cumhazard == "log"){
                ## a * exp +/-b = exp(log(a) +/- b)
                ## formula 4.10 p 58 (Beyersmann et al. 2012)
                object$cumhazard.lower <- object$cumhazard * exp(-zval * object$cumhazard.se)
                object$cumhazard.upper <- object$cumhazard * exp(+zval * object$cumhazard.se)
            }
        }

        if("survival" %in% type){
            if(object$transform.survival == "none"){
                object$survival.lower <- object$survival.upper <- matrix(NA, nrow = NROW(object$survival), ncol = NCOL(object$survival)) 
                object$survival.lower[] <- apply(object$survival - zval * object$survival.se,2,pmax,0) # to keep matrix format even when object$survival contains only one line
                object$survival.upper[] <- apply(object$survival + zval * object$survival.se,2,pmin,1) # to keep matrix format even when object$survival contains only one line
            }else if(object$transform.survival == "log"){
                ## a * exp +/-b = exp(log(a) +/- b)
                ## formula 4.10 p 58 (Beyersmann et al. 2012)
                object$survival.lower <- object$survival.upper <- matrix(NA, nrow = NROW(object$survival), ncol = NCOL(object$survival)) 
                object$survival.lower[] <- apply(object$survival * exp(-zval * object$survival.se),2,pmax,0) # to keep matrix format even when object$survival contains only one line
                object$survival.upper[] <- object$survival * exp(+zval * object$survival.se) # to keep matrix format even when object$survival contains only one line
            }else if(object$transform.survival == "loglog"){
                ## exp(-exp(log(-log(a)) +/- b)) = exp(-exp(log(-log(a)))exp(+/- b)) = exp(-(-log(a))exp(+/- b)) = exp(log(a)exp(+/- b)) = a ^ exp(+/- b)
                ## formula 4.16 p 59 (Beyersmann et al. 2012)
                object$survival.lower <- object$survival^(exp( + zval * object$survival.se ))
                object$survival.upper <- object$survival^(exp( - zval * object$survival.se ))                
            }else if(object$transform.survival == "cloglog"){
                ## formula 4.21 p 62 (Beyersmann et al. 2012)
                object$survival.lower <- 1 - (1-object$survival)^(exp( - zval * object$survival.se ))
                object$survival.upper <- 1 - (1-object$survival)^(exp( + zval * object$survival.se ))                

            }
        }
        
    }

    ## confidence bands
    if(object$band && nsim.band > 0){
        if ("cumhazard" %in% type){
            if(!is.na(seed)){set.seed(seed)}
            object$quantile.band <- confBandCox(iid = object$cumhazard.iid,
                                                se = object$cumhazard.se,
                                                n.sim = nsim.band,
                                                conf.level = conf.level)
            
            seBand_tempo <- colMultiply_cpp(object$cumhazard.se,object$quantile.band)

            if("cumhazard" %in% type){
                if(object$transform.cumhazard == "none"){
                    object$cumhazard.lowerBand <- matrix(NA, nrow = NROW(object$cumhazard), ncol = NCOL(object$cumhazard)) # to keep matrix format even when object$cumhazard contains only one line
                    object$cumhazard.lowerBand[] <- apply(object$cumhazard - zval * seBand_tempo,2,pmax,0) # to keep matrix format even when object$cumhazard contains only one line
                    object$cumhazard.upperBand <- object$cumhazard + zval * seBand_tempo
                }else if(object$transform.cumhazard == "log"){
                    ## a * exp +/-b = exp(log(a) +/- b)
                    ## formula 4.10 p 58 (Beyersmann et al. 2012)
                    object$cumhazard.lowerBand <- object$cumhazard * exp(-zval * seBand_tempo)
                    object$cumhazard.upperBand <- object$cumhazard * exp(+zval * seBand_tempo)
                }
            } 
        }

        if("survival" %in% type){
            if(!is.na(seed)){set.seed(seed)}
            object$quantile.band <- confBandCox(iid = object$survival.iid,
                                                se = object$survival.se,
                                                n.sim = nsim.band,
                                                conf.level = conf.level)
            
            seBand_tempo <- colMultiply_cpp(object$survival.se,object$quantile.band)/zval

            if(object$transform.survival == "none"){
                object$survival.lowerBand <- object$survival.upperBand <- matrix(NA, nrow = NROW(object$survival), ncol = NCOL(object$survival)) 
                object$survival.lowerBand[] <- apply(object$survival - zval * seBand_tempo,2,pmax,0) # to keep matrix format even when object$survival contains only one line
                object$survival.upperBand[] <- apply(object$survival + zval * seBand_tempo,2,pmin,1) # to keep matrix format even when object$survival contains only one line
            }else if(object$transform.survival == "log"){
                ## a * exp +/-b = exp(log(a) +/- b)
                ## formula 4.10 p 58 (Beyersmann et al. 2012)
                object$survival.lowerBand <- object$survival.upperBand <- matrix(NA, nrow = NROW(object$survival), ncol = NCOL(object$survival)) 
                object$survival.lowerBand[] <- apply(object$survival * exp(-zval * seBand_tempo),2,pmax,0) # to keep matrix format even when object$survival contains only one line
                object$survival.upperBand[] <- object$survival * exp(+zval * seBand_tempo) # to keep matrix format even when object$survival contains only one line
            }else if(object$transform.survival == "loglog"){
                ## exp(-exp(log(-log(a)) +/- b)) = exp(-exp(log(-log(a)))exp(+/- b)) = exp(-(-log(a))exp(+/- b)) = exp(log(a)exp(+/- b)) = a ^ exp(+/- b)
                ## formula 4.16 p 59 (Beyersmann et al. 2012)
                object$survival.lowerBand <- object$survival^(exp( + zval * seBand_tempo ))
                object$survival.upperBand <- object$survival^(exp( - zval * seBand_tempo ))                
            }else if(object$transform.survival == "cloglog"){
                ## formula 4.21 p 62 (Beyersmann et al. 2012)
                object$survival.lowerBand <- 1 - (1-object$survival)^(exp( - zval * seBand_tempo ))
                object$survival.upperBand <- 1 - (1-object$survival)^(exp( + zval * seBand_tempo ))                

            }
        }

    }

    ## check NA
    indexNA <- union(which(is.na(object[[paste0(type,".se")]])),
                     which(is.nan(object[[paste0(type,".se")]])))
    if(length(indexNA)>0){

        if(se){
            object[[paste0(type,".lower")]][indexNA] <- NA
            object[[paste0(type,".upper")]][indexNA] <- NA
        }
        if(band){
            indexNA2 <- union(which(rowSums(is.na(object[[paste0(type,".se")]]))>0),
                              which(rowSums(is.nan(object[[paste0(type,".se")]]))>0))
            object[["quantile.band"]][indexNA2] <- NA
            object[[paste0(type,".lowerBand")]][indexNA2,] <- NA
            object[[paste0(type,".upperBand")]][indexNA2,] <- NA
        }
        
    }

    
    ## export
    return(object)
}

######################################################################
### confint.predictCox.R ends here
