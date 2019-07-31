### confint.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: jul 31 2019 (11:15) 
##           By: Brice Ozenne
##     Update #: 553
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
##' @param nsim.band [integer, >0]the number of simulations used to compute the quantiles for the confidence bands.
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
##' @param bootci.method [character] Method for constructing bootstrap confidence intervals.
##' Either "perc" (the default), "norm", "basic", "stud", or "bca".
##' @param ... not used.
##'
##' @details
##' Confidence bands and confidence intervals computed via the influence function
##' are automatically restricted to the interval [0;1]. \cr \cr
##'
##' Confidence intervals obtained via bootstrap are computed
##' using the \code{boot.ci} function of the \code{boot} package.
##' p-value are obtained using test inversion method
##' (finding the smallest confidence level such that the interval contain the null hypothesis).
##' 
##' @author Brice Ozenne

## * confint.ate (examples)
##' @rdname confint.ate
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
##' print(fit.ate, type = "meanRisk")
##' 
##' 
##' ## manual calculation of se
##' dd <- copy(d)
##' dd$X1 <- rep(factor("T0", levels = paste0("T",0:2)), NROW(dd))
##' out <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, average.iid = TRUE)
##' term1 <- -out$survival.average.iid
##' term2 <- sweep(1-out$survival, MARGIN = 2, FUN = "-", STATS = colMeans(1-out$survival))
##' sqrt(colSums((term1 + term2/NROW(d))^2)) 
##' ## fit.ate$meanRisk[Treatment=="T0",meanRisk.se]
##' 
##' ## note
##' out2 <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, iid = TRUE)
##' mean(out2$survival.iid[,1,1])
##' out$survival.average.iid[1,1]
##' 
##' ## check confidence intervals (no transformation)
##' fit.ate$meanRisk[, .(lower = meanRisk + qnorm(0.025) * meanRisk.se,
##'                      upper = meanRisk + qnorm(0.975) * meanRisk.se)]
##' ## fit.ate$meanRisk[,cbind(meanRisk.lower,meanRisk.upper)]
##' 
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' outCI <- confint(fit.ate,
##'                  meanRisk.transform = "loglog", diffRisk.transform = "atanh",
##'                  ratioRisk.transform = "log")
##' print(outCI, type = "meanRisk")
##' 
##' newse <- fit.ate$meanRisk[, meanRisk.se/(meanRisk*log(meanRisk))]
##' fit.ate$meanRisk[, .(lower = exp(-exp(log(-log(meanRisk)) - 1.96 * newse)),
##'                      upper = exp(-exp(log(-log(meanRisk)) + 1.96 * newse)))]
##' ## outCI$meanRisk[,cbind(meanRisk.lower,meanRisk.upper)]
## * confint.ate (code)
##' @rdname confint.ate
##' @method confint ate
##' @export
confint.ate <- function(object,
                        parm = NULL,
                        level = 0.95,
                        nsim.band = 1e4,
                        meanRisk.transform = "none",
                        diffRisk.transform = "none",
                        ratioRisk.transform = "none",
                        seed = NA,
                        bootci.method = "perc",
                        ...){

    if(object$se[[1]] == FALSE && object$band[[1]] == FALSE){
        message("No confidence interval is computed \n",
                "Set argument \'se\' to TRUE when calling ate \n")
        return(object)
    }

    ## ** hard copy
    ## needed otherwise meanRisk and riskComparison are modified in the original object
    object$meanRisk <- data.table::copy(object$meanRisk)
    object$riskComparison <- data.table::copy(object$riskComparison)

    ## ** compute CI
    if(!is.null(object$boot)){
        object <- confintBoot.ate(object,
                                  bootci.method = bootci.method,
                                  conf.level = level,
                                  ...)

    }else{
        object <- confintIID.ate(object,
                                 nsim.band = nsim.band,
                                 meanRisk.transform = meanRisk.transform,
                                 diffRisk.transform = diffRisk.transform,
                                 ratioRisk.transform = ratioRisk.transform,
                                 seed = seed,
                                 conf.level = level,
                                 ...)
    }
    object$meanRisk[] ## ensure direct print
    object$riskComparison[] ## ensure direct print

    ## ** export
    class(object) <- "ate"
    return(object)
}

## * confintBoot.ate
confintBoot.ate <- function(object,
                            bootci.method,
                            conf.level){

    valid.boot <- c("norm","basic","stud","perc","wald","quantile")
    bootci.method <- match.arg(bootci.method, valid.boot)
    ## normalize arguments
    bootci.method <- tolower(bootci.method) ## convert to lower case
    name.estimate <- names(object$boot$t0)
    n.estimate <- length(name.estimate)
    index <- 1:n.estimate
    alpha <- 1-conf.level
    
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
    
    ## store arguments in the object
    object$conf.level <- conf.level
    object$bootci.method <- bootci.method

    ## real number of bootstrap samples
    test.NA <- !is.na(object$boot$t)
    test.Inf <- !is.infinite(object$boot$t)
    n.boot <- colSums(test.NA*test.Inf)

    ## standard error
    boot.se <- sqrt(apply(object$boot$t, 2, var, na.rm = TRUE))
    boot.mean <- colMeans(object$boot$t, na.rm = TRUE)

    ## confidence interval
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
            out <- boot::boot.ci(object$boot,
                                 conf = conf.level,
                                 type = bootci.method,
                                 index = iP)[[slot.boot.ci]][index.lowerCI:index.upperCI]
            return(setNames(out,c("lower","upper")))
        }    
    }),silent=TRUE)
    if (class(try.CI)[1]=="try-error"){
        warning("Could not construct bootstrap confidence limits")
        boot.CI <- matrix(rep(NA,2*length(index)),ncol=2)
    } else{
        boot.CI <- do.call(rbind,ls.CI)
    }
    
    ## pvalue
    null <- setNames(rep(NA,length(name.estimate)),name.estimate)
    null[grep("^compRisk:diff", name.estimate)] <- 0
    null[grep("^compRisk:ratio", name.estimate)] <- 1

    boot.p <- sapply(index, function(iP){ # iP <- 1
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
                                       boot::boot.ci(object$boot,
                                                     conf = 1-p.value,
                                                     type = bootci.method,
                                                     index = iP)[[slot.boot.ci]][side.CI]
                                   })
            return(p.value)
        }
    })

    ## group
    dt.tempo <- data.table(name = name.estimate, mean = boot.mean, se = boot.se, boot.CI, p.value = boot.p)
    keep.col <- setdiff(names(dt.tempo),"name")
    
    vecNames.meanRisk <- paste0("meanRisk",c(".bootstrap",".se",".lower",".upper"))
    object$meanRisk[,c(vecNames.meanRisk) := dt.tempo[grep("^meanRisk",dt.tempo$name),.SD,.SDcols = setdiff(keep.col,"p.value")]]
    
    vecNames.diffRisk <- paste0("diff",c(".bootstrap",".se",".lower",".upper",".p.value"))
    object$riskComparison[,c(vecNames.diffRisk) := dt.tempo[grep("^compRisk:diff",dt.tempo$name),.SD,.SDcols = keep.col]]

    vecNames.ratioRisk <- paste0("ratio",c(".bootstrap",".se",".lower",".upper",".p.value"))
    object$riskComparison[,c(vecNames.ratioRisk) := dt.tempo[grep("^compRisk:ratio",dt.tempo$name),.SD,.SDcols = keep.col]]
    
    return(object)    
}

## * confintIID.ate 
confintIID.ate <- function(object,
                           conf.level,
                           nsim.band,
                           meanRisk.transform,
                           diffRisk.transform,
                           ratioRisk.transform,
                           seed){

    if(object$se[[1]] == FALSE && object$band[[1]] == FALSE){
        message("No confidence interval/band computed \n",
                "Set argument \'se\' or argument \'band\' to TRUE when calling predictCSC \n")
        return(object)
    }

    ## ** check arguments
    if(is.null(object$iid)){
        stop("Cannot re-compute standard error or confidence bands without the iid decomposition \n",
             "Set argument \'iid\' to TRUE when calling ate \n")
    }
    object$meanRisk.transform <- match.arg(meanRisk.transform, c("none","log","loglog","cloglog"))
    object$diffRisk.transform <- match.arg(diffRisk.transform, c("none","atanh"))
    object$ratioRisk.transform <- match.arg(ratioRisk.transform, c("none","log"))

    ## ** prepare
    times <- object$times
    n.times <- length(times)
    contrasts <- object$contrasts
    n.contrasts <- length(contrasts)
    allContrasts <- utils::combn(contrasts, m = 2)
    n.allContrasts <- NCOL(allContrasts)
    n.obs <- NROW(object$iid[[1]])
    
    ## ** meanRisk
    
    ## initialize
    if(object$se){
        object$meanRisk[,
                        c("meanRisk.se","meanRisk.lower","meanRisk.upper") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }
    if(object$band){
        object$meanRisk[,
                        c("meanRisk.quantileBand","meanRisk.lowerBand","meanRisk.upperBand") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }

    ## reshape data
    estimate.mR <- matrix(NA, nrow = n.contrasts, ncol = n.times)
    iid.mR <- array(NA, dim = c(n.contrasts, n.times, n.obs))
    for(iC in 1:n.contrasts){ ## iT <- 1
        estimate.mR[iC,] <- object$meanRisk[object$meanRisk[[1]]==contrasts[iC],
                                            .SD$meanRisk]
        iid.mR[iC,,] <- t(object$iid[[contrasts[iC]]])
    }
    se.mR <- sqrt(apply(iid.mR^2, MARGIN = 1:2, sum))

    ## compute
    CIBP.mR <- transformCIBP(estimate = estimate.mR,
                             se = se.mR,
                             iid = iid.mR,
                             null = NA,
                             conf.level = conf.level,
                             nsim.band = nsim.band,
                             seed = seed,
                             type = object$meanRisk.transform,
                             min.value = switch(object$meanRisk.transform,
                                                "none" = 0,
                                                "log" = NULL,
                                                "loglog" = NULL,
                                                "cloglog" = NULL),
                             max.value = switch(object$meanRisk.transform,
                                                "none" = 1,
                                                "log" = 1,
                                                "loglog" = NULL,
                                                "cloglog" = NULL),
                             ci = object$se,
                             band = object$band,
                             p.value = FALSE)

    ## store
    for(iC in 1:n.contrasts){ ## iT <- 1
        if(object$se){
            object$meanRisk[object$meanRisk[[1]] == contrasts[iC],
                            c("meanRisk.se","meanRisk.lower","meanRisk.upper") := list(se.mR[iC,],CIBP.mR$lower[iC,],CIBP.mR$upper[iC,])]
        }
        if(object$band){
            object$meanRisk[object$meanRisk[[1]] == contrasts[iC],
                            c("meanRisk.quantileBand","meanRisk.lowerBand","meanRisk.upperBand") := list(CIBP.mR$quantileBand[iC],CIBP.mR$lowerBand[iC,],CIBP.mR$upperBand[iC,])]
        }
    }
    
    ## ** diffRisk: se, CI/CB
    ## initialize
    if(object$se){
        object$riskComparison[,
                              c("diff.se","diff.lower","diff.upper","diff.p.value") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }
    if(object$band){
        object$riskComparison[,
                              c("diff.quantileBand","diff.lowerBand","diff.upperBand") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }

    ## reshape data
    estimate.dR <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    iid.dR <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
    
    for(iC in 1:n.allContrasts){ ## iC <- 1
        estimate.dR[iC,] <- object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                                  .SD$diff]
        iid.dR[iC,,] <- t(object$iid[[allContrasts[2,iC]]] - object$iid[[allContrasts[1,iC]]])
    }
    se.dR <- sqrt(apply(iid.dR^2, MARGIN = 1:2, sum))

    ## compute
    CIBP.dR <- transformCIBP(estimate = estimate.dR,
                             se = se.dR,
                             iid = iid.dR,
                             null = 0,
                             conf.level = conf.level,
                             nsim.band = nsim.band,
                             seed = seed,
                             type = object$diffRisk.transform,
                             min.value = switch(object$diffRisk.transform,
                                                "none" = -1,
                                                "atanh" = NULL),
                             max.value = switch(object$diffRisk.transform,
                                                "none" = 1,
                                                "atanh" = NULL),
                             ci = object$se,
                             band = object$band,
                             p.value = object$se)


    ## store
    for(iC in 1:n.allContrasts){ ## iT <- 1
        if(object$se){
            object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                  c("diff.se","diff.lower","diff.upper","diff.p.value") := list(se.dR[iC,],CIBP.dR$lower[iC,],CIBP.dR$upper[iC,],CIBP.dR$p.value[iC,])]
        }
        if(object$band){
            object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                  c("diff.quantileBand","diff.lowerBand","diff.upperBand") := list(CIBP.dR$quantileBand[iC],CIBP.dR$lowerBand[iC,],CIBP.dR$upperBand[iC,])]
        }
    }

    ## ** ratioRisk: se, CI/CB
    ## initialize
    if(object$se){
        object$riskComparison[,
                              c("ratio.se","ratio.lower","ratio.upper","ratio.p.value") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }
    if(object$band){
        object$riskComparison[,
                              c("ratio.quantileBand","ratio.lowerBand","ratio.upperBand") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA))]
    }

    ## reshape data
    estimate.rR <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    iid.rR <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
    
    for(iC in 1:n.allContrasts){ ## iC <- 1
        estimate.rR[iC,] <- object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                                  .SD$ratio]
        term1 <- sweep(object$iid[[allContrasts[2,iC]]], MARGIN = 2, FUN = "/",
                       STATS = object$meanRisk[object$meanRisk[[1]] == allContrasts[1,iC],.SD$meanRisk])
        term2 <- sweep(object$iid[[allContrasts[1,iC]]], MARGIN = 2, FUN = "*",
                       STATS = object$meanRisk[object$meanRisk[[1]] == allContrasts[2,iC],.SD$meanRisk]/object$meanRisk[object$meanRisk[[1]] == allContrasts[1,iC],.SD$meanRisk]^2)
        iid.rR[iC,,] <- t(term1 - term2)
    }
    se.rR <- sqrt(apply(iid.rR^2, MARGIN = 1:2, sum))
    
    ## compute
    CIBP.rR <- transformCIBP(estimate = estimate.rR,
                             se = se.rR,
                             iid = iid.rR,
                             null = 1,
                             conf.level = conf.level,
                             nsim.band = nsim.band,
                             seed = seed,
                             type = object$ratioRisk.transform,
                             min.value = switch(object$ratioRisk.transform,
                                                "none" = 0,
                                                "log" = NULL),
                             max.value = NULL,
                             ci = object$se,
                             band = object$band,
                             p.value = object$se)

    ## store
    for(iC in 1:n.allContrasts){ ## iT <- 1
        if(object$se){
            object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                  c("ratio.se","ratio.lower","ratio.upper","ratio.p.value") := list(se.rR[iC,],CIBP.rR$lower[iC,],CIBP.rR$upper[iC,],CIBP.rR$p.value[iC,])]
        }
        if(object$band){
            object$riskComparison[(object$riskComparison[[1]] == allContrasts[1,iC]) & (object$riskComparison[[2]] == allContrasts[2,iC]),
                                  c("ratio.quantileBand","ratio.lowerBand","ratio.upperBand") := list(CIBP.rR$quantileBand[iC],CIBP.rR$lowerBand[iC,],CIBP.rR$upperBand[iC,])]
        }
    }

    ## ** re-order columns
    suffix <- c("",".se")
    if(object$se){
        suffix <- c(suffix, c(".p.value",".lower",".upper"))
    }
    if(object$band){
        suffix <- c(suffix, c(".quantileBand",".lowerBand",".upperBand"))
    }
    data.table::setcolorder(object$riskComparison,
                            neworder = c(names(object$riskComparison)[1:3],
                                         paste0("diff",suffix),
                                         paste0("ratio",suffix)
                                         )
                            )
    ## ** export
    object$conf.level <- conf.level
    object$nsim.band <- nsim.band
    return(object)
}

######################################################################
### confint.ate.R ends here
