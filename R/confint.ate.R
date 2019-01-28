### confint.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: Jan 28 2019 (16:33) 
##           By: Thomas Alexander Gerds
##     Update #: 478
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
##' ## manual calculation of se
##' dd <- copy(d)
##' dd$X1 <- rep(factor("T0", levels = paste0("T",0:2)), NROW(dd))
##' out <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, average.iid = TRUE)
##' term1 <- -out$survival.average.iid
##' term2 <- sweep(1-out$survival, MARGIN = 2, FUN = "-", STATS = colMeans(1-out$survival))
##' sqrt(colSums((term1 + term2/NROW(d))^2)) 
##' 
##' ## note
##' out2 <- predictCox(fit, newdata = dd, se = TRUE, times = 1:3, iid = TRUE)
##' mean(out2$survival.iid[,1,1])
##' out$survival.average.iid[1,1]
##' 
##' ## check confidence intervals (no transformation)
##' fit.ate$meanRisk[, .(lower = meanRisk + qnorm(0.025) * meanRisk.se,
##'                      upper = meanRisk + qnorm(0.975) * meanRisk.se)]
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

    if(object$se == FALSE && object$band == FALSE){
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
    ls.CI <- lapply(index, function(iP){ # iP <- 1
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
    })
    boot.CI <- do.call(rbind,ls.CI)
    
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


    if(object$se == FALSE && object$band == FALSE){
        message("No confidence interval/band computed \n",
                "Set argument \'se\' or argument \'band\' to TRUE when calling predictCSC \n")
        return(object)
    }

    ## ** check arguments
    if(object$band && (is.null(object$meanRisk$meanRisk.se) || is.null(object$riskComparison$diff.se) || is.null(object$riskComparison$ratio.se)) ){
        stop("Cannot compute confidence bands \n",
             "Set argument \'se\' to TRUE when calling ate \n")
    }
    if(object$band && (is.null(object$meanRisk.iid) || is.null(object$diffRisk.iid) || is.null(object$ratioRisk.iid))){
        stop("Cannot compute confidence bands \n",
             "Set argument \'iid\' to TRUE when calling ate \n")
    }
    object$meanRisk.transform <- match.arg(meanRisk.transform, c("none","log","loglog","cloglog"))
    object$diffRisk.transform <- match.arg(diffRisk.transform, c("none","atanh"))
    object$ratioRisk.transform <- match.arg(ratioRisk.transform, c("none","log"))

    ## ** prepare
    n.times <- length(object$times)
    n.treatment <- length(object$contrasts)
    allContrasts <- utils::combn(object$contrasts, m = 2)
    n.allContrasts <- NCOL(allContrasts)
    n.obs <- NCOL(object$meanRisk.iid)

    ## ** meanRisk: se, CI/CB

    ## reshape data
    ls.index.tempo <- list()
    estimate.tempo <- matrix(NA, nrow = n.treatment, ncol = n.times)
    se.tempo <- matrix(NA, nrow = n.treatment, ncol = n.times)
    if(object$band){
        iid.tempo <- array(NA, dim = c(n.treatment, n.times, n.obs))
    }else{
        iid.tempo <- NULL
    }
    
    for(iT in 1:n.treatment){ ## iT <- 1
        ls.index.tempo[[iT]] <- which(object$meanRisk[[1]]==object$contrasts[iT])

        estimate.tempo[iT,] <- object$meanRisk[["meanRisk"]][ls.index.tempo[[iT]]]
        se.tempo[iT,] <- object$meanRisk[["meanRisk.se"]][ls.index.tempo[[iT]]]
        if(object$band){
            iid.tempo[iT,,] <- object$meanRisk.iid[ls.index.tempo[[iT]],]
        }        
    }

    ## compute
    outCIBP.meanRisk <- transformCIBP(estimate = estimate.tempo,
                                      se = se.tempo,
                                      iid = iid.tempo,
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
    vec.index.tempo <- unlist(ls.index.tempo)
    if(object$se){
        object$meanRisk <- cbind(object$meanRisk[vec.index.tempo],
                                 meanRisk.lower=as.double(t(outCIBP.meanRisk$lower)),
                                 meanRisk.upper=as.double(t(outCIBP.meanRisk$upper)))
    }
    if(object$band){
        object$meanRisk <- cbind(object$meanRisk[vec.index.tempo],
                                 meanRisk.quantileBand=rep(as.double(t(outCIBP.meanRisk$quantileBand)),length.out=length(vec.index.tempo)), 
                                 meanRisk.lowerBand=as.double(t(outCIBP.meanRisk$lowerBand)),
                                 meanRisk.upperBand=as.double(t(outCIBP.meanRisk$upperBand)))
    }
    ## ** diffRisk: se, CI/CB

    ## reshape data
    ls.index.tempo <- list()
    estimate.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    se.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    if(object$band){
        iid.tempo <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
    }else{
        iid.tempo <- NULL
    }
    
    for(iC in 1:n.allContrasts){ ## iC <- 1
        ls.index.tempo[[iC]] <- intersect(which(object$riskComparison[[1]]==allContrasts[1,iC]),
                                          which(object$riskComparison[[2]]==allContrasts[2,iC]))

        estimate.tempo[iC,] <- object$riskComparison[["diff"]][ls.index.tempo[[iC]]]
        se.tempo[iC,] <- object$riskComparison[["diff.se"]][ls.index.tempo[[iC]]]
        if(object$band){
            iid.tempo[iC,,] <- object$diffRisk.iid[ls.index.tempo[[iC]],]
        }        
    }

    ## compute
    outCIBP.diffRisk <- transformCIBP(estimate = estimate.tempo,
                                      se = se.tempo,
                                      iid = iid.tempo,
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
    vec.index.tempo <- unlist(ls.index.tempo)
    if(object$se){
        object$riskComparison[,c("diff.lower","diff.upper","diff.p.value")] <- as.numeric(NA)
        object$riskComparison[vec.index.tempo, c("diff.lower") := as.double(t(outCIBP.diffRisk$lower))]
        object$riskComparison[vec.index.tempo, c("diff.upper") := as.double(t(outCIBP.diffRisk$upper))]
        object$riskComparison[vec.index.tempo, c("diff.p.value") := as.double(t(outCIBP.diffRisk$p.value))]
    }
    if(object$band){
        object$riskComparison[,c("diff.quantileBand","diff.lowerBand","diff.upperBand")] <- as.numeric(NA)
        object$riskComparison[vec.index.tempo, c("diff.quantileBand") := rep(as.double(t(outCIBP.diffRisk$quantileBand)),length.out=length(vec.index.tempo))]
        object$riskComparison[vec.index.tempo, c("diff.lowerBand") := as.double(t(outCIBP.diffRisk$lowerBand))]
        object$riskComparison[vec.index.tempo, c("diff.upperBand") := as.double(t(outCIBP.diffRisk$upperBand))]
    }
    ## ** ratioRisk: se, CI/CB

    ## reshape data
    ls.index.tempo <- list()
    estimate.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    se.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
    if(object$band){
        iid.tempo <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
    }else{
        iid.tempo <- NULL
    }
    
    for(iC in 1:n.allContrasts){ ## iC <- 1
        ls.index.tempo[[iC]] <- intersect(which(object$riskComparison[[1]]==allContrasts[1,iC]),
                            which(object$riskComparison[[2]]==allContrasts[2,iC]))

        estimate.tempo[iC,] <- object$riskComparison[["ratio"]][ls.index.tempo[[iC]]]
        se.tempo[iC,] <- object$riskComparison[["ratio.se"]][ls.index.tempo[[iC]]]
        if(object$band){
            iid.tempo[iC,,] <- object$ratioRisk.iid[ls.index.tempo[[iC]],]
        }        
    }

    ## compute
    outCIBP.ratioRisk <- transformCIBP(estimate = estimate.tempo,
                                       se = se.tempo,
                                       iid = iid.tempo,
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
    vec.index.tempo <- unlist(ls.index.tempo)
    if(object$se){
        object$riskComparison[,c("ratio.lower","ratio.upper","ratio.p.value")] <- as.numeric(NA)
        object$riskComparison[vec.index.tempo, c("ratio.lower") := as.double(t(outCIBP.ratioRisk$lower))]
        object$riskComparison[vec.index.tempo, c("ratio.upper") := as.double(t(outCIBP.ratioRisk$upper))]
        object$riskComparison[vec.index.tempo, c("ratio.p.value") := as.double(t(outCIBP.ratioRisk$p.value))]
    }
    if(object$band){
        object$riskComparison[,c("ratio.quantileBand","ratio.lowerBand","ratio.upperBand")] <- as.numeric(NA)
        object$riskComparison[vec.index.tempo, c("ratio.quantileBand") := rep(as.double(t(outCIBP.ratioRisk$quantileBand)),length.out=length(vec.index.tempo))]
        object$riskComparison[vec.index.tempo, c("ratio.lowerBand") := as.double(t(outCIBP.ratioRisk$lowerBand))]
        object$riskComparison[vec.index.tempo, c("ratio.upperBand") := as.double(t(outCIBP.ratioRisk$upperBand))]

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
