### confint.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: Mar  7 2022 (08:28) 
##           By: Thomas Alexander Gerds
##     Update #: 1005
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.ate (documentation)
##' @title Confidence Intervals and Confidence Bands for the Predicted Absolute Risk (Cumulative Incidence Function)
##' @description Confidence intervals and confidence Bands for the predicted absolute risk (cumulative incidence function).
##' @name confint.ate
##' 
##' @param object A \code{ate} object, i.e. output of the \code{ate} function.
##' @param parm not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] Level of confidence.
##' @param n.sim [integer, >0] the number of simulations used to compute the quantiles for the confidence bands and/or perform adjustment for multiple comparisons.
##' @param estimator [character] The type of estimator relative to which the estimates should be displayed. 
##' @param contrasts [character vector] levels of the treatment variable for which the risks should be assessed and compared. Default is to consider all levels.
##' @param allContrasts [2-row character matrix] levels of the treatment variable to be compared. Default is to consider all pairwise comparisons.
##' @param meanRisk.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the mean risk in small samples.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}.
##' @param diffRisk.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the risk difference in small samples.
##' Can be \code{"none"}, \code{"atanh"}.
##' @param ratioRisk.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the risk ratio in small samples.
##' Can be \code{"none"}, \code{"log"}.
##' @param seed [integer, >0] seed number set when performing simulation for the confidence bands.
##' If not given or NA no seed is set.
##' @param ci [logical] should the confidence intervals be computed?
##' @param band [logical] should the confidence bands be computed?
##' @param p.value [logical] should the p-values/adjusted p-values be computed?
##' Requires argument \code{ci} and/or \code{band} to be \code{TRUE}.
##' @param method.band [character] method used to adjust for multiple comparisons.
##' Can be any element of \code{p.adjust.methods} (e.g. \code{"holm"}), \code{"maxT-integration"}, or \code{"maxT-simulation"}. 
##' @param alternative [character] a character string specifying the alternative hypothesis,
##' must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}.
##' @param bootci.method [character] Method for constructing bootstrap confidence intervals.
##' Either "perc" (the default), "norm", "basic", "stud", or "bca".
##' @param ... not used.
##'
##' @details
##' Argument \code{ci}, \code{band}, \code{p.value}, \code{method.band}, \code{alternative}, \code{meanRisk.transform}, \code{diffRisk.transform}, \code{ratioRisk.transform} are only active when the \code{ate} object contains the influence function.
##' Argument \code{bootci.method} is only active when the \code{ate} object contains bootstrap samples.
##' 
##' \strong{Influence function}: confidence bands and confidence intervals computed via the influence function are automatically restricted to the interval of definition of the parameter (e.g. [0;1] for the average risk).
##' Single step max adjustment for multiple comparisons, i.e. accounting for the correlation between the test statistics but not for the ordering of the tests, can be performed setting the arguemnt \code{method.band} to \code{"maxT-integration"} or \code{"maxT-simulation"}. The former uses numerical integration (\code{pmvnorm} and \code{qmvnorm} to perform the adjustment while the latter using simulation. Both assume that the test statistics are jointly normally distributed.
##'
##' \strong{Bootstrap}: confidence intervals obtained via bootstrap are computed
##' using the \code{boot.ci} function of the \code{boot} package.
##' p-value are obtained using test inversion method
##' (finding the smallest confidence level such that the interval contain the null hypothesis).
##' 
##' @author Brice Ozenne

## * confint.ate (examples)
##' @examples
##' library(survival)
##' library(data.table)
##' 
##' ## ## generate data ####
##' set.seed(10)
##' d <- sampleData(70,outcome="survival")
##' d[, X1 := paste0("T",rbinom(.N, size = 2, prob = c(0.51)))]
##' ## table(d$X1)
##' 
##' #### stratified Cox model ####
##' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
##'              data=d, ties="breslow", x = TRUE, y = TRUE)
##' 
##' #### average treatment effect ####
##' fit.ate <- ate(fit, treatment = "X1", times = 1:3, data = d,
##'                se = TRUE, iid = TRUE, band = TRUE)
##' summary(fit.ate)
##' dt.ate <- as.data.table(fit.ate)
##' 
##' ## manual calculation of se
##' dd <- copy(d)
##' dd$X1 <- rep(factor("T0", levels = paste0("T",0:2)), NROW(dd))
##' out <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, average.iid = TRUE)
##' term1 <- -out$survival.average.iid
##' term2 <- sweep(1-out$survival, MARGIN = 2, FUN = "-", STATS = colMeans(1-out$survival))
##' sqrt(colSums((term1 + term2/NROW(d))^2)) 
##' ## fit.ate$meanRisk[treatment=="T0",se]
##' 
##' ## note
##' out2 <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, iid = TRUE)
##' mean(out2$survival.iid[1,1,])
##' out$survival.average.iid[1,1]
##' 
##' ## check confidence intervals (no transformation)
##' dt.ate[,.(lower = pmax(0,estimate + qnorm(0.025) * se),
##'           lower2 = lower,
##'           upper = estimate + qnorm(0.975) * se,
##'           upper2 = upper)]
##' 
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' outCI <- confint(fit.ate,
##'                  meanRisk.transform = "loglog", diffRisk.transform = "atanh",
##'                  ratioRisk.transform = "log")
##' summary(outCI, type = "risk", short = TRUE)
##' 
##' dt.ate[type == "meanRisk", newse := se/(estimate*log(estimate))]
##' dt.ate[type == "meanRisk", .(lower = exp(-exp(log(-log(estimate)) - 1.96 * newse)),
##'                         upper = exp(-exp(log(-log(estimate)) + 1.96 * newse)))]

## * confint.ate (code)
##' @rdname confint.ate
##' @method confint ate
##' @export
confint.ate <- function(object,
                        parm = NULL,
                        level = 0.95,
                        n.sim = 1e4,
                        estimator = object$estimator,
                        contrasts = object$contrasts,
                        allContrasts = object$allContrasts,
                        meanRisk.transform = "none",
                        diffRisk.transform = "none",
                        ratioRisk.transform = "none",
                        seed = NA,
                        ci = object$inference$se,
                        band = object$inference$band,
                        p.value = TRUE,
                        method.band = "maxT-simulation",
                        alternative = "two.sided",
                        bootci.method = "perc",
                        ...){

    ## ** check arguments
    if(object$inference$se == FALSE && object$inference$band == FALSE){
        message("No confidence interval is computed \n",
                "Set argument \'se\' to TRUE when calling the ate function \n")
        return(object)
    }
    if(any(contrasts %in% object$contrasts == FALSE)){
        stop("Incorrect values for the argument \'contrasts\' \n",
             "Possible values: \"",paste(object$contrasts,collapse="\" \""),"\" \n")
    }
    if(!is.null(allContrasts)){
        if(any(allContrasts %in% object$contrasts == FALSE)){
            stop("Incorrect values for the argument \'allContrasts\' \n",
                 "Possible values: \"",paste(object$contrasts,collapse="\" \""),"\" \n")
        }
        if(any(interaction(allContrasts[1,],allContrasts[2,]) %in% interaction(object$allContrasts[1,],object$allContrasts[2,]) == FALSE)){
            stop("Incorrect combinations for the argument \'allContrasts\' \n",
                 "Possible combinations: \"",paste(interaction(object$allContrasts[1,],object$allContrasts[2,]),collapse="\" \""),"\" \n")
        }
        if(!is.matrix(allContrasts) || NROW(allContrasts)!=2){
            stop("Argument \'allContrasts\' must be a matrix with 2 rows \n")
        }
    }
    estimator <- match.arg(estimator, object$estimator, several.ok = TRUE)
    
    meanRisk.transform <- match.arg(meanRisk.transform, c("none","log","loglog","cloglog"))
    diffRisk.transform <- match.arg(diffRisk.transform, c("none","atanh"))
    ratioRisk.transform <- match.arg(ratioRisk.transform, c("none","log"))

    bootci.method <- match.arg(bootci.method, c("norm","basic","stud","perc","wald","quantile"))
    bootstrap <- object$inference$bootstrap
    
    ## ** initialize
    rm.col <- c("se","lower","upper","p.value","quantileBand","lowerBand","upperBand","adj.p.value")

    out <- list(meanRisk = data.table::copy(object$meanRisk[,.SD,.SDcols = setdiff(names(object$meanRisk),rm.col)]),
                diffRisk = data.table::copy(object$diffRisk[,.SD,.SDcols = setdiff(names(object$diffRisk),rm.col)]),
                ratioRisk = data.table::copy(object$ratioRisk[,.SD,.SDcols = setdiff(names(object$ratioRisk),rm.col)]),
                inference = object$inference,
                inference.allContrasts = allContrasts,
                inference.contrasts = contrasts,
                transform = c("meanRisk" = meanRisk.transform,
                              "diffRisk" = diffRisk.transform,
                              "ratioRisk" = ratioRisk.transform))
    out$inference["conf.level"] <- level
    out$inference["ci"] <- ci
    out$inference["p.value"] <- p.value
    out$inference["alternative"] <- alternative
    out$inference["band"] <- band
    if(bootstrap){
        out$inference["bootci.method"] <- bootci.method
    }else if(band){
        out$inference["method.band"] <- method.band
        out$inference["n.sim"] <- n.sim
    }
    
    out$meanRisk[,c("se") := as.numeric(NA)]
    out$diffRisk[,c("se") := as.numeric(NA)]
    out$ratioRisk[,c("se") := as.numeric(NA)]

    if(ci){
        out$meanRisk[,c("lower","upper") := as.numeric(NA)]
        out$diffRisk[,c("lower","upper",if(p.value){"p.value"}) := as.numeric(NA)]
        out$ratioRisk[,c("lower","upper",if(p.value){"p.value"}) := as.numeric(NA)]
    }
    if(band>0){
        if(method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            out$meanRisk[,c("quantileBand","lowerBand","upperBand") := as.numeric(NA)]
            out$diffRisk[,c("quantileBand","lowerBand","upperBand") := as.numeric(NA)]
            out$ratioRisk[,c("quantileBand","lowerBand","upperBand") := as.numeric(NA)]
        }
        if(p.value){
            out$diffRisk[,c("adj.p.value") := as.numeric(NA)]
            out$ratioRisk[,c("adj.p.value") := as.numeric(NA)]
        }
    }

    ## ** compute CI
    if(bootstrap){
        if(band && !identical(method.band,"none")){
            stop("Adjustment for multiple comparisons not implemented for the boostrap. \n")
        }
        if(!identical(alternative,"two.sided")){
            stop("Only two sided tests are implemented for the boostrap. \n")
        }
        if(!identical(out$inference.contrasts,object$contrasts)){
            stop("Argument \'contrast\' is not available when using boostrap. \n")
        }
        if(!identical(out$inference.allContrasts,object$allContrasts)){
            stop("Argument \'allContrast\' is not available when using boostrap. \n")
        }
        out <- confintBoot.ate(object, estimator = estimator, out = out, seed = seed)
    }else{
        out <- confintIID.ate(object, estimator = estimator, out = out, seed = seed)
    }
    ## out$meanRisk[] ## ensure direct print
    ## out$diffRisk[] ## ensure direct print
    ## out$ratioRisk[] ## ensure direct print

    ## ** export
    return(out)
}

## * confintBoot.ate
confintBoot.ate <- function(object, estimator, out, seed){

    name.estimate <- names(object$boot$t0)
    n.estimate <- length(name.estimate)
    contrasts <- out$inference.contrasts
    n.contrasts <- length(contrasts)
    allContrasts <- out$inference.allContrasts
    n.allContrasts <- NCOL(allContrasts)
    collapse.allContrasts <- interaction(allContrasts[1,],allContrasts[2,])
    index <- 1:n.estimate
    conf.level <- out$inference$conf.level
    alpha <- 1-conf.level
    bootci.method <- out$inference$bootci.method
    p.value <- out$inference$p.value
    ci <- out$inference$ci
    slot.boot.ci <- switch(bootci.method,
                           "norm" = "normal",
                           "basic" = "basic",
                           "stud" = "student",
                           "perc" = "percent",
                           "bca" = "bca")
    index.lowerCI <- switch(bootci.method,
                            "norm" = 2,
                            "basic" = 4,
                            "stud" = 4,
                            "perc" = 4,
                            "bca" = 4)
    index.upperCI <- switch(bootci.method,
                            "norm" = 3,
                            "basic" = 5,
                            "stud" = 5,
                            "perc" = 5,
                            "bca" = 5)

    ## real number of bootstrap samples
    test.NA <- !is.na(object$boot$t)
    test.Inf <- !is.infinite(object$boot$t)
    n.boot <- colSums(test.NA*test.Inf)

    ## standard error
    boot.se <- sqrt(apply(object$boot$t, 2, var, na.rm = TRUE))
    boot.mean <- colMeans(object$boot$t, na.rm = TRUE)
    
    ## confidence interval
    if(ci){
        try.CI <- try(ls.CI <- lapply(index, function(iP){ # iP <- 1        
            if(n.boot[iP]==0){
                return(c(lower = NA, upper = NA))
            }else if(bootci.method == "wald"){
                return(c(lower = as.double(object$boot$t0[iP] + qnorm(alpha/2) * boot.se[iP]),
                         upper = as.double(object$boot$t0[iP] - qnorm(alpha/2) * boot.se[iP])
                         ))
            }else if(bootci.method == "quantile"){
                return(c(lower = as.double(quantile(object$boot$t[,iP], probs = alpha/2, na.rm = TRUE)),
                         upper = as.double(quantile(object$boot$t[,iP], probs = 1-(alpha/2), na.rm = TRUE))
                         ))
            }else{
                if (requireNamespace("boot",quietly=TRUE)){
                    out <- boot::boot.ci(object$boot,
                                         conf = conf.level,
                                         type = bootci.method,
                                         index = iP)[[slot.boot.ci]][index.lowerCI:index.upperCI]
                }else{
                    stop("Package 'boot' requested to obtain confidence intervals, but not installed.")
                }
                return(setNames(out,c("lower","upper")))
            }    
        }),silent=TRUE)
        if (inherits(x=(try.CI),what="try-error")){
            warning("Could not construct bootstrap confidence limits")
            boot.CI <- matrix(rep(NA,2*length(index)),ncol=2)
        } else{
            boot.CI <- do.call(rbind,ls.CI)
        }
    }

    ## pvalue
    if(p.value){
    null <- setNames(rep(NA,length(name.estimate)),name.estimate)
    null[grep("^diff", name.estimate)] <- 0
    null[grep("^ratio", name.estimate)] <- 1

    boot.p <- sapply(index, function(iP){ # iP <- 25
        iNull <- null[iP]
        if(is.na(iNull)){return(NA)}
        iEstimate <- object$boot$t0[iP]
        iSE <- boot.se[iP]
            
        if(n.boot[iP]==0){
            return(NA)
        }else if(bootci.method == "wald"){
            return(2*(1-stats::pnorm(abs((iEstimate-iNull)/iSE))))
        }else if(bootci.method == "quantile"){
            if(iEstimate > iNull){
                return(mean(object$boot$t[,iP] > iNull))
            }else if(iEstimate < iNull){
                return(mean(object$boot$t[,iP] < iNull))
            }else{
                return(1)
            }
        }else{
            ## search confidence level such that quantile of CI which is close to 0
            p.value <- boot2pvalue(x = object$boot$t[,iP],
                                   null = iNull,
                                   estimate = iEstimate,
                                   alternative = "two.sided",
                                   FUN.ci = function(p.value, sign.estimate, ...){ ## p.value <- 0.4
                                       side.CI <- c(index.lowerCI,index.upperCI)[2-sign.estimate]
                                       if (requireNamespace("boot",quietly=TRUE)){
                                           boot::boot.ci(object$boot,
                                                         conf = 1-p.value,
                                                         type = bootci.method,
                                                         index = iP)[[slot.boot.ci]][side.CI]
                                       }else{
                                           stop("Package 'boot' requested to obtain confidence intervals, but not installed.")
                                       }
                                   })
            return(p.value)
        }
    })
    }
    
    ## store
    vcov.meanRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)
    vcov.diffRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)
    vcov.ratioRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)

    for(iE in 1:length(estimator)){ ## iE <- 1
        iEstimator <- estimator[iE]

        indexMean <- grep(paste0("^mean.",iEstimator),name.estimate)
        indexDiff <- grep(paste0("^difference.",iEstimator),name.estimate)
        indexRatio <- grep(paste0("^ratio.",iEstimator),name.estimate)

        out$meanRisk[iEstimator, c("estimate.boot") :=  boot.mean[indexMean], on = "estimator"]
        out$diffRisk[iEstimator, c("estimate.boot") :=  boot.mean[indexDiff], on = "estimator"]
        out$ratioRisk[iEstimator, c("estimate.boot") :=  boot.mean[indexRatio], on = "estimator"]

        out$meanRisk[iEstimator, c("se") :=  boot.se[indexMean], on = "estimator"]
        out$diffRisk[iEstimator, c("se") :=  boot.se[indexDiff], on = "estimator"]
        out$ratioRisk[iEstimator, c("se") :=  boot.se[indexRatio], on = "estimator"]
        
        vcov.meanRisk[[iEstimator]] <- setNames(lapply(contrasts, function(iC){ ## iC <- contrasts[1]
            var(object$boot$t[,grep(paste0(iC,"$"),name.estimate[indexMean], value = TRUE),drop=FALSE])
        }), contrasts)
        vcov.diffRisk[[iEstimator]] <- setNames(lapply(collapse.allContrasts, function(iC){ ## iC <- collapse.allContrasts[1]
            var(object$boot$t[,grep(paste0(iC,"$"),name.estimate[indexDiff], value = TRUE),drop=FALSE])
        }), collapse.allContrasts)
        vcov.ratioRisk[[iEstimator]] <- setNames(lapply(collapse.allContrasts, function(iC){ ## iC <- collapse.allContrasts[1]
            var(object$boot$t[,grep(paste0(iC,"$"),name.estimate[indexRatio], value = TRUE),drop=FALSE])
        }), collapse.allContrasts)
        
        
        if(ci){
            out$meanRisk[iEstimator, c("lower","upper") :=  list(boot.CI[indexMean,"lower"],boot.CI[indexMean,"upper"]), on = "estimator"]
            out$diffRisk[iEstimator, c("lower","upper") :=  list(boot.CI[indexDiff,"lower"],boot.CI[indexDiff,"upper"]), on = "estimator"]
            out$ratioRisk[iEstimator, c("lower","upper") :=  list(boot.CI[indexRatio,"lower"],boot.CI[indexRatio,"upper"]), on = "estimator"]
        }

        if(p.value){
            out$diffRisk[iEstimator, c("p.value") :=  boot.p[indexDiff], on = "estimator"]
            out$ratioRisk[iEstimator, c("p.value") :=  boot.p[indexRatio], on = "estimator"]
        }
    }
    ## ** export
    if(attr(object$estimator,"TD")){
        setcolorder(out$meanRisk, neworder = c("estimator","time","landmark","treatment","estimate","estimate.boot","se","lower","upper"))
        setcolorder(out$diffRisk, neworder = c("estimator","time","landmark","A","B","estimate.A","estimate.B","estimate","estimate.boot","se","lower","upper","p.value"))
        setcolorder(out$ratioRisk, neworder = c("estimator","time","landmark","A","B","estimate.A","estimate.B","estimate","estimate.boot","se","lower","upper","p.value"))
    }else{
        setcolorder(out$meanRisk, neworder = c("estimator","time","treatment","estimate","estimate.boot","se","lower","upper"))
        setcolorder(out$diffRisk, neworder = c("estimator","time","A","B","estimate.A","estimate.B","estimate","estimate.boot","se","lower","upper","p.value"))
        setcolorder(out$ratioRisk, neworder = c("estimator","time","A","B","estimate.A","estimate.B","estimate","estimate.boot","se","lower","upper","p.value"))
    }
    data.table::setattr(out$meanRisk, name = "vcov", value = vcov.meanRisk)
    data.table::setattr(out$diffRisk, name = "vcov", value = vcov.diffRisk)
    data.table::setattr(out$ratioRisk, name = "vcov", value = vcov.ratioRisk)
    return(out)    
}

## * confintIID.ate 
confintIID.ate <- function(object, estimator, out, seed){
    n.estimator <- length(estimator)

    ## ** check arguments
    if(is.null(object$iid) || any(sapply(object$iid[estimator],is.null))){
        stop("Cannot re-compute standard error or confidence bands without the iid decomposition \n",
             "Set argument \'iid\' to TRUE when calling the ate function \n")
    }
                             
    ## ** prepare
    times <- object$eval.times
    n.times <- length(times)
    contrasts <- out$inference.contrasts
    n.contrasts <- length(contrasts)
    allContrasts <- out$inference.allContrasts
    n.allContrasts <- NCOL(allContrasts)
    collapse.allContrasts <- interaction(allContrasts[1,],allContrasts[2,])
    n.obs <- max(stats::na.omit(object$n[-1]))
    conf.level <- out$inference$conf.level
    alternative <- out$inference$alternative
    ci <- out$inference$ci
    band <- out$inference$band
    method.band <- out$inference$method.band
    p.value <- out$inference$p.value
    n.sim <- out$inference$n.sim

    meanRisk.transform <- out$transform["meanRisk"]
    diffRisk.transform <- out$transform["diffRisk"]
    ratioRisk.transform <- out$transform["ratioRisk"]

    vcov.meanRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)
    vcov.diffRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)
    vcov.ratioRisk <- setNames(vector(mode = "list", length = length(estimator)), estimator)
    
    ## ** run
    for(iE in 1:n.estimator){ ## iE <- 1
        iEstimator <- estimator[iE]
        ## ** meanRisk
        ## cat("meanRisk \n")

        ## reshape data
        estimate.mR <- matrix(NA, nrow = n.contrasts, ncol = n.times)
        iid.mR <- array(NA, dim = c(n.obs, n.times, n.contrasts))
        for(iC in 1:n.contrasts){ ## iC <- 1
            estimate.mR[iC,] <- object$meanRisk[list(iEstimator,contrasts[iC]), .SD$estimate, on = c("estimator","treatment")]
            iid.mR[,,iC] <- object$iid[[iEstimator]][[contrasts[iC]]]
        }
        se.mR <- t(sqrt(apply(iid.mR^2, MARGIN = 2:3, sum)))

        ## compute
        CIBP.mR <- transformCIBP(estimate = estimate.mR,
                                 se = se.mR,
                                 iid = iid.mR,
                                 null = NA,
                                 conf.level = conf.level,
                                 alternative = "two.sided",
                                 n.sim = n.sim,
                                 seed = seed,
                                 type = meanRisk.transform,
                                 min.value = switch(meanRisk.transform,
                                                    "none" = 0,
                                                    "log" = NULL,
                                                    "loglog" = NULL,
                                                    "cloglog" = NULL),
                                 max.value = switch(meanRisk.transform,
                                                    "none" = 1,
                                                    "log" = 1,
                                                    "loglog" = NULL,
                                                    "cloglog" = NULL),
                                 ci = ci,
                                 band = band,
                                 method.band = method.band,
                                 p.value = FALSE)

        ## store
        vcov.meanRisk[[iEstimator]] <- setNames(vector(mode = "list", length = n.contrasts), contrasts)
        for(iC in 1:n.contrasts){ ## iC <- 2
            iRowIndex <- which((out$meanRisk$estimator==iEstimator)*(out$meanRisk$treatment==contrasts[iC])==1)
            
            out$meanRisk[iRowIndex, c("se") := se.mR[iC,]]
            if(ci){
                out$meanRisk[iRowIndex, c("lower","upper") := list(CIBP.mR$lower[iC,],CIBP.mR$upper[iC,])]
            }
            if((band>0) && (method.band %in% c("bonferroni","maxT-integration","maxT-simulation"))){
                out$meanRisk[iRowIndex, c("quantileBand","lowerBand","upperBand") := list(CIBP.mR$quantileBand[iC],CIBP.mR$lowerBand[iC,],CIBP.mR$upperBand[iC,])]
            }
            if(n.times==1){
                vcov.meanRisk[[iEstimator]][[iC]] <- sum(iid.mR[,,iC]^2)
            }else{
                vcov.meanRisk[[iEstimator]][[iC]] <- crossprod(iid.mR[,,iC])
            }
        }

        
        

        ## ** diffRisk: se, CI/CB
        ## cat("diffRisk \n")

        ## reshape data
        estimate.dR <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
        iid.dR <- array(NA, dim = c(n.obs, n.times, n.allContrasts))
        for(iC in 1:n.allContrasts){ ## iC <- 1
            estimate.dR[iC,] <- object$diffRisk[list(iEstimator,allContrasts[1,iC],allContrasts[2,iC]), .SD$estimate, on = c("estimator","A","B")]
            iid.dR[,,iC] <- object$iid[[iEstimator]][[allContrasts[2,iC]]] - object$iid[[iEstimator]][[allContrasts[1,iC]]]
        }
        se.dR <- t(sqrt(apply(iid.dR^2, MARGIN = 2:3, sum)))

        ## compute
        CIBP.dR <- transformCIBP(estimate = estimate.dR,
                                 se = se.dR,
                                 iid = iid.dR,
                                 null = 0,
                                 conf.level = conf.level,
                                 alternative = alternative,
                                 n.sim = n.sim,
                                 seed = seed,
                                 type = diffRisk.transform,
                                 min.value = switch(diffRisk.transform,
                                                    "none" = -1,
                                                    "atanh" = NULL),
                                 max.value = switch(diffRisk.transform,
                                                    "none" = 1,
                                                    "atanh" = NULL),
                                 ci = ci,
                                 band = band,
                                 method.band = method.band,
                                 p.value = p.value)
        
        ## store
        vcov.diffRisk[[iEstimator]] <- setNames(vector(mode = "list", length = n.allContrasts), collapse.allContrasts)
        for(iC in 1:n.allContrasts){ ## iC <- 1
            iRowIndex <- which((out$diffRisk$estimator==iEstimator)*(out$diffRisk$A==allContrasts[1,iC])*(out$diffRisk$B==allContrasts[2,iC])==1)
            out$diffRisk[iRowIndex, c("se") := se.dR[iC,]]
            if(ci){
                out$diffRisk[iRowIndex, c("lower","upper") := list(CIBP.dR$lower[iC,],CIBP.dR$upper[iC,])]
                if(p.value){
                    out$diffRisk[iRowIndex, c("p.value") := CIBP.dR$p.value[iC,]]
                }
            }
            if(band>0){
                if(method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
                    out$diffRisk[iRowIndex, c("quantileBand","lowerBand","upperBand") := list(CIBP.dR$quantileBand[iC],CIBP.dR$lowerBand[iC,],CIBP.dR$upperBand[iC,])]
                }
                if(p.value){
                    out$diffRisk[iRowIndex, c("adj.p.value") := CIBP.dR$adj.p.value[iC,]]
                }
            }
            if(n.times==1){
                vcov.diffRisk[[iEstimator]][[iC]] <- sum(iid.dR[,,iC]^2)
            }else{
                vcov.diffRisk[[iEstimator]][[iC]] <- crossprod(iid.dR[,,iC])
            }
        }

        ## ** ratioRisk: se, CI/CB
        ## cat("ratioRisk \n")

        ## reshape data
        estimate.rR <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
        iid.rR <- array(NA, dim = c(n.obs, n.times, n.allContrasts))
    
        for(iC in 1:n.allContrasts){ ## iC <- 1
            estimate.rR[iC,] <- object$ratioRisk[list(iEstimator,allContrasts[1,iC],allContrasts[2,iC]), .SD$estimate, on = c("estimator","A","B")]

            factor1 <- object$meanRisk[list(iEstimator, allContrasts[1,iC]),.SD$estimate, on = c("estimator","treatment")]
            factor2 <- object$meanRisk[list(iEstimator, allContrasts[2,iC]),.SD$estimate, on = c("estimator","treatment")]

            term1 <- rowMultiply_cpp(object$iid[[iEstimator]][[allContrasts[2,iC]]],
                                     scale = 1/factor1)
            term2 <- rowMultiply_cpp(object$iid[[iEstimator]][[allContrasts[1,iC]]],
                                     scale = factor2/factor1^2)
            iid.rR[,,iC] <- term1 - term2
        }
        se.rR <- t(sqrt(apply(iid.rR^2, MARGIN = 2:3, sum)))

        ## compute
        CIBP.rR <- transformCIBP(estimate = estimate.rR,
                                 se = se.rR,
                                 iid = iid.rR,
                                 null = 1,
                                 conf.level = conf.level,
                                 alternative = alternative,
                                 n.sim = n.sim,
                                 seed = seed,
                                 type = ratioRisk.transform,
                                 min.value = switch(ratioRisk.transform,
                                                    "none" = 0,
                                                    "log" = NULL),
                                 max.value = NULL,
                                 ci = ci,
                                 band = band,
                                 method.band = method.band,
                                 p.value = p.value)

        ## store
        vcov.ratioRisk[[iEstimator]] <- setNames(vector(mode = "list", length = n.allContrasts), collapse.allContrasts)
        for(iC in 1:n.allContrasts){ ## iC <- 1
            iRowIndex <- which((out$ratioRisk$estimator==iEstimator)*(out$ratioRisk$A==allContrasts[1,iC])*(out$ratioRisk$B==allContrasts[2,iC])==1)

            out$ratioRisk[iRowIndex, c("se") := se.rR[iC,]]
            if(ci){
                out$ratioRisk[iRowIndex, c("lower","upper") := list(CIBP.rR$lower[iC,],CIBP.rR$upper[iC,])]
                if(p.value){
                    out$ratioRisk[iRowIndex, c("p.value") := CIBP.rR$p.value[iC,]]
                }
            }
            if(band>0){
                if(method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
                    out$ratioRisk[iRowIndex, c("quantileBand","lowerBand","upperBand") := list(CIBP.rR$quantileBand[iC],CIBP.rR$lowerBand[iC,],CIBP.rR$upperBand[iC,])]
                }
                if(p.value){
                    out$ratioRisk[iRowIndex, c("adj.p.value") := CIBP.rR$adj.p.value[iC,]]
                }
            }
            if(n.times==1){
                vcov.ratioRisk[[iEstimator]][[iC]] <- sum(iid.rR[,,iC]^2)
            }else{
                vcov.ratioRisk[[iEstimator]][[iC]] <- crossprod(iid.rR[,,iC])
            }
        }
    }

    ## ** export
    data.table::setattr(out$meanRisk, name = "vcov", value = vcov.meanRisk)
    data.table::setattr(out$diffRisk, name = "vcov", value = vcov.diffRisk)
    data.table::setattr(out$ratioRisk, name = "vcov", value = vcov.ratioRisk)
    return(out)
}

######################################################################
### confint.ate.R ends here
