### confint.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (18:06) 
##           By: Brice Ozenne
##     Update #: 328
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
##' @param type.boot [character] Method for constructing bootstrap confidence intervals.
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
##' d <- sampleData(40,outcome="survival")
##' 
##' #### stratified Cox model ####
##' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
##'              data=d, ties="breslow", x = TRUE, y = TRUE)
##'
##' #### average treatment effect ####
##' fit.pred <- ate(fit, treatment = "X1", times = 1:3, data = d)
##'
##' ## manual calculation of se
##' dd <- copy(d)
##' dd$X1 <- factor(0, levels = 0:1)
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
##' ## add confidence intervals computed on the original scale
##' confint(fit.pred,
##' meanRisk.transform = "none", diffRisk.transform = "none", ratioRisk.transform = "none"
##' )
##' 
##' fit.pred$meanRisk[, .(lower = meanRisk - 1.96 * se,
##'                       upper = meanRisk + 1.96 * se)]
##'
##' ## add confidence intervals computed on the log-log scale
##' ## and backtransformed
##' confint(fit.pred,
##' meanRisk.transform = "loglog", diffRisk.transform = "atanh", ratioRisk.transform = "log"
##' )
##' 
##' newse <- fit.pred$meanRisk[, se/(meanRisk*log(meanRisk))]
##' fit.pred$meanRisk[, .(lower = exp(-exp(log(-log(meanRisk)) - 1.96 * newse)),
##'                       upper = exp(-exp(log(-log(meanRisk)) + 1.96 * newse)))]

## * confint.ate (code)
##' @rdname confint.ate
##' @export
confint.ate <- function(object,
                        parm = NULL,
                        level = 0.95,
                        nsim.band = 1e4,
                        meanRisk.transform = "none",
                        diffRisk.transform = "none",
                        ratioRisk.transform = "none",
                        seed = NA,
                        type.boot = "perc",
                        ...){


    ## compute confidence intervals and p.values in matrix form
    if(object$se || object$band){
        ## hard copy
        object$meanRisk <- data.table::copy(object$meanRisk)
        object$riskComparison <- data.table::copy(object$riskComparison)

        if(!is.null(object$boot)){
            object <- confintBoot.ate(object,
                                      type.boot = type.boot,
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
    }else{
        message("No confidence interval is computed \n",
                "Set argument \'se\' to TRUE when calling ate \n")
    }

    ## export
    return(object)
}

## * confintBoot.ate
confintBoot.ate <- function(object,
                            type.boot,
                            conf.level){

    ## normalize arguments
    type.boot <- tolower(type.boot) ## convert to lower case
    name.estimate <- names(object$boot$t0)
    n.estimate <- length(name.estimate)
    index <- 1:n.estimate
    alpha <- 1-conf.level
    
    slot.boot.ci <- switch(type.boot,
                           "norm" = "normal",
                           "basic" = "basic",
                           "stud" = "student",
                           "perc" = "percent",
                           "bca" = "bca")
    index.lowerCI <- switch(type.boot,
                            "norm" = 2,
                            "basic" = 4,
                            "stud" = 4,
                            "perc" = 4,
                            "bca" = 4)
    index.upperCI <- switch(type.boot,
                            "norm" = 3,
                            "basic" = 5,
                            "stud" = 5,
                            "perc" = 5,
                            "bca" = 5)
    
    ## store arguments in the object
    object$conf.level <- conf.level
    object$type.boot <- type.boot

    ## real number of bootstrap samples
    test.NA <- !is.na(object$boot$t)
    test.Inf <- !is.infinite(object$boot$t)
    n.boot <- colSums(test.NA*test.Inf)

    ## standard error
    boot.se <- sqrt(apply(object$boot$t, 2, var, na.rm = TRUE))

    ## confidence interval
    ls.CI <- lapply(index, function(iP){ # iP <- 1
        if(n.boot[iP]==0){
            return(c(lower = NA, upper = NA))
        }else if(type.boot == "wald"){
            return(c(lower = as.double(object$boot$t0[iP] + qnorm(alpha/2) * boot.se[iP]),
                     upper = as.double(object$boot$t0[iP] - qnorm(alpha/2) * boot.se[iP])
                     ))
        }else if(type.boot == "quantile"){
            return(c(lower = as.double(quantile(object$boot$t[,iP], probs = alpha/2, na.rm = TRUE)),
                     upper = as.double(quantile(object$boot$t[,iP], probs = 1-(alpha/2), na.rm = TRUE))
                     ))
        }else{
            out <- boot::boot.ci(object$boot,
                                 conf = conf.level,
                                 type = type.boot,
                                 index = iP)[[slot.boot.ci]][index.lowerCI:index.upperCI]
            return(setNames(out,c("lower","upper")))
        }    
    })
    boot.CI <- do.call(rbind,ls.CI)
    
    ## pvalue
    null <- setNames(rep(0,length(name.estimate)),name.estimate)
    null[grep("^compRisk:ratio", name.estimate)] <- 1

    boot.p <- sapply(index, function(iP){ # iP <- 1
        iNull <- null[iP]
        iEstimate <- object$boot$t0[iP]
        iSE <- boot.se[iP]
            
        if(n.boot[iP]==0){
            return(NA)
        }else if(type.boot == "wald"){
            return(2*(1-stats::pnorm(abs((iEstimate-iNull)/iSE))))
        }else if(type.boot == "quantile"){
            if(iEstimate>iNull){
                return(mean(object$boot$t[,iP] > iNull))
            }else if(iEstimate<null){
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
                                                     type = type.boot,
                                                     index = iP)[[slot.boot.ci]][side.CI]
                                   })
            return(p.value)
        }
    })

    ## group
    object$confint <- data.table(name.estimate, boot.p, boot.se, boot.CI)
    return(object)    
}

## * confintIId.ate 
confintIID.ate <- function(object,
                           conf.level,
                           nsim.band,
                           meanRisk.transform,
                           diffRisk.transform,
                           ratioRisk.transform,
                           seed){

    .I <- NULL ## [:CRANcheck:] data.table
    
    ## ** check arguments
    if(!is.null(object$meanRisk.transform) && object$meanRisk.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    if(!is.null(object$diffRisk.transform) && object$diffRisk.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    if(!is.null(object$diffRatio.transform) && object$diffRatio.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    object$meanRisk.transform <- match.arg(meanRisk.transform, c("none","log","loglog","cloglog"))
    object$diffRisk.transform <- match.arg(diffRisk.transform, c("none","atanh"))
    object$ratioRisk.transform <- match.arg(ratioRisk.transform, c("none","log"))

    if(object$band){
        n.times <- length(object$times)
        n.treatment <- length(object$contrasts)
        allContrasts <- utils::combn(object$contrasts, m = 2)
        n.allContrasts <- NCOL(allContrasts)
        name.allContrasts <- apply(allContrasts,2,paste, collapse = ".")
        name.repAllContrasts <- paste0(object$riskComparison[[1]],".",object$riskComparison[[2]])
        n.obs <- NCOL(object$ate.iid)
    }

    ## ** quantile
    zval <- stats::qnorm(1 - (1-conf.level)/2, mean = 0, sd = 1)

    ## ** ate
    if(object$meanRisk.transform != "none"){
        ## transform standard error
        object$meanRisk[,c("se") := transformSE(estimate = object$meanRisk[["meanRisk"]],
                                                se = object$meanRisk[["se"]],
                                                type = object$meanRisk.transform)]

        if(object$band){
            object$ate.iid <- transformIID(estimate = object$meanRisk[["meanRisk"]],
                                           iid = object$ate.iid,
                                           type = object$meanRisk.transform,
                                           format = "matrix")
        }

    }
    
    ate.min.value <- switch(object$meanRisk.transform,
                            "none" = 0,
                            "log" = NULL,
                            "loglog" = NULL,
                            "cloglog" = NULL)
    ate.max.value <- switch(object$meanRisk.transform,
                            "none" = 1,
                            "log" = 1,
                            "loglog" = NULL,
                            "cloglog" = NULL)
    ## confidence intervals
    if(object$se){
        object$meanRisk[,c("lower","upper") := transformCI(estimate = object$meanRisk[["meanRisk"]],
                                                           se = object$meanRisk[["se"]],
                                                           quantile = zval,
                                                           type = object$meanRisk.transform,
                                                           format = "vector",
                                                           min.value = ate.min.value,
                                                           max.value = ate.max.value)]

        indexNA <- object$meanRisk[,.I[is.na(.SD$meanRisk)+is.na(.SD$se)+is.nan(.SD$meanRisk)+is.nan(.SD$se)>0]]
        if(length(indexNA)>0){
            object$meanRisk[indexNA, c("lower","upper") := NA]
        }
    }

    ## confidence bands
    if(object$band){
        ## find quantiles for the bands
        if(!is.na(seed)){set.seed(seed)}

        iid.tempo <- array(NA, dim = c(n.treatment, n.times, n.obs))
        se.tempo <- matrix(NA, nrow = n.treatment, ncol = n.times)
        for(iT in 1:n.treatment){ ## iT <- 1
            iIndex <- which(object$meanRisk[[1]]==object$contrasts[iT])
            iid.tempo[iT,,] <- object$ate.iid[iIndex,]
            se.tempo[iT,] <- object$meanRisk[["se"]][iIndex]
        }
        
        quantileBand.tempo <- confBandCox(iid = iid.tempo,
                                          se = se.tempo,
                                          n.sim = nsim.band,
                                          conf.level = conf.level)
        quantileBand.tempo <- setNames(quantileBand.tempo,object$contrast)
            
        object$meanRisk[, c("quantileBand") := quantileBand.tempo[.SD[[1]]]]
        object$meanRisk[,c("lowerBand","upperBand") := transformCI(estimate = object$meanRisk[["meanRisk"]],
                                                                   se = object$meanRisk[["se"]],
                                                                   quantile = object$meanRisk[["quantileBand"]],
                                                                   type = object$meanRisk.transform,
                                                                   format = "vector",
                                                                   min.value = ate.min.value,
                                                                   max.value = ate.max.value)]

        indexNA <- object$meanRisk[,is.na(.SD$meanRisk)+is.na(.SD$se)+is.nan(.SD$meanRisk)+is.nan(.SD$se)>0]
        indexTNA <- which(tapply(indexNA,object$meanRisk[[1]],sum)>0)
        if(length(indexTNA)>0){            
            object$meanRisk[which(object$meanRisk[[1]] %in% object$contrasts[indexTNA]),
                            c("lowerBand","upperBand") := NA]
        }
    }
 
    ## ** diff treatment
    if(object$diffRisk.transform != "none"){
        ## transform standard error
        object$riskComparison[,c("diff.se") := transformSE(estimate = object$riskComparison[["diff"]],
                                                           se = object$riskComparison[["diff.se"]],
                                                           type = object$diffRisk.transform)]

    }
    diffAte.min.value <- switch(object$diffRisk.transform,
                                "none" = -1,
                                "atanh" = NULL)
    diffAte.max.value <- switch(object$diffRisk.transform,
                                "none" = 1,
                                "atanh" = NULL)
    ## confidence intervals
    if(object$se){
        object$riskComparison[,c("diff.lower","diff.upper") := transformCI(estimate = object$riskComparison[["diff"]],
                                                                           se = object$riskComparison[["diff.se"]],
                                                                           quantile = zval,
                                                                           type = object$diffRisk.transform,
                                                                           format = "vector",
                                                                           min.value = diffAte.min.value,
                                                                           max.value = diffAte.max.value)]
        ## p.value
        object$riskComparison[,c("diff.p.value") := transformP(estimate = object$riskComparison[["diff"]],
                                                               se = object$riskComparison[["diff.se"]],
                                                               null = 0,
                                                               type = object$diffRisk.transform)]

        ## NA
        indexNA <- object$riskComparison[,.I[is.na(.SD$diff)+is.na(.SD$diff.se)+is.nan(.SD$diff)+is.nan(.SD$diff.se)>0]]
        if(length(indexNA)>0){
            object$riskComparison[indexNA, c("diff.lower","diff.upper","diff.p.value") := NA]
        }
    }

    ## confidence bands
    if(object$band){
        ## find quantiles for the bands
        iid.tempo <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
        se.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
        for(iC in 1:n.allContrasts){ ## iC <- 1
            iIndex <- intersect(which(object$riskComparison[[1]]==allContrasts[1,iC]),
                                which(object$riskComparison[[2]]==allContrasts[2,iC]))
            iid.tempo[iC,,] <- object$diffAte.iid[iIndex,]
            se.tempo[iC,] <- object$riskComparison[["diff.se"]][iIndex]
        }
        
        quantileBand.tempo <- confBandCox(iid = iid.tempo,
                                          se = se.tempo,
                                          n.sim = nsim.band,
                                          conf.level = conf.level)
        quantileBand.tempo <- setNames(quantileBand.tempo,name.allContrasts)
        
        object$riskComparison[, c("diff.quantileBand") := quantileBand.tempo[paste0(.SD[[1]],".",.SD[[2]])]]
        object$riskComparison[,c("diff.lowerBand","diff.upperBand") := transformCI(estimate = object$riskComparison[["diff"]],
                                                                                   se = object$riskComparison[["diff.se"]],
                                                                                   quantile = object$riskComparison[["diff.quantileBand"]],
                                                                                   type = object$diffRisk.transform,
                                                                                   format = "vector",
                                                                                   min.value = diffAte.min.value,
                                                                                   max.value = diffAte.max.value)]

        indexNA <- object$riskComparison[,is.na(.SD$diff)+is.na(.SD$diff.se)+is.nan(.SD$diff)+is.nan(.SD$diff.se)>0]
        indexTNA <- which(tapply(indexNA,name.repAllContrasts,sum)>0)
        if(length(indexTNA)>0){            
            object$riskComparison[which(name.repAllContrasts %in% name.allContrasts[indexTNA]),
                                  c("diff.lowerBand","diff.upperBand") := NA]
        }
    }
    
    ## ** ratio treatment
    if(object$ratioRisk.transform != "none"){
        ## transform standard error
        object$riskComparison[,c("ratio.se") := transformSE(estimate = object$riskComparison[["ratio"]],
                                                           se = object$riskComparison[["ratio.se"]],
                                                           type = object$ratioRisk.transform)]

    }
    ratioAte.min.value <- switch(object$ratioRisk.transform,
                                 "none" = 0,
                                 "log" = NULL)
    ratioAte.max.value <- NULL
    
    ## confidence intervals
    if(object$se){
        object$riskComparison[,c("ratio.lower","ratio.upper") := transformCI(estimate = object$riskComparison[["ratio"]],
                                                                             se = object$riskComparison[["ratio.se"]],
                                                                             quantile = zval,
                                                                             type = object$ratioRisk.transform,
                                                                             format = "vector",
                                                                             min.value = ratioAte.min.value,
                                                                             max.value = ratioAte.max.value)]
        ## p.value
        object$riskComparison[,c("ratio.p.value") := transformP(estimate = object$riskComparison[["ratio"]],
                                                                se = object$riskComparison[["ratio.se"]],
                                                                null = 1,
                                                                type = object$ratioRisk.transform)]

        ## NA
        indexNA <- object$riskComparison[,.I[is.na(.SD$ratio)+is.na(.SD$ratio.se)+is.nan(.SD$ratio)+is.nan(.SD$ratio.se)>0]]
        if(length(indexNA)>0){
            object$riskComparison[indexNA, c("ratio.lower","ratio.upper","ratio.p.value") := NA]
        }
    }

    if(object$band){
        ## find quantiles for the bands
        iid.tempo <- array(NA, dim = c(n.allContrasts, n.times, n.obs))
        se.tempo <- matrix(NA, nrow = n.allContrasts, ncol = n.times)
        for(iC in 1:n.allContrasts){ ## iC <- 1
            iIndex <- intersect(which(object$riskComparison[[1]]==allContrasts[1,iC]),
                                which(object$riskComparison[[2]]==allContrasts[2,iC]))
            iid.tempo[iC,,] <- object$ratioAte.iid[iIndex,]
            se.tempo[iC,] <- object$riskComparison[["ratio.se"]][iIndex]
        }
        
        quantileBand.tempo <- confBandCox(iid = iid.tempo,
                                          se = se.tempo,
                                          n.sim = nsim.band,
                                          conf.level = conf.level)
        quantileBand.tempo <- setNames(quantileBand.tempo,name.allContrasts)
        
        object$riskComparison[, c("ratio.quantileBand") := quantileBand.tempo[paste0(.SD[[1]],".",.SD[[2]])]]
        object$riskComparison[,c("ratio.lowerBand","ratio.upperBand") := transformCI(estimate = object$riskComparison[["ratio"]],
                                                                                     se = object$riskComparison[["ratio.se"]],
                                                                                     quantile = object$riskComparison[["ratio.quantileBand"]],
                                                                                     type = object$ratioRisk.transform,
                                                                                     format = "vector",
                                                                                     min.value = ratioAte.min.value,
                                                                                     max.value = ratioAte.max.value)]

        indexNA <- object$riskComparison[,is.na(.SD$ratio)+is.na(.SD$ratio.se)+is.nan(.SD$ratio)+is.nan(.SD$ratio.se)>0]
        indexTNA <- which(tapply(indexNA,name.repAllContrasts,sum)>0)
        if(length(indexTNA)>0){            
            object$riskComparison[which(name.repAllContrasts %in% name.allContrasts[indexTNA]),
                                  c("ratio.lowerBand","ratio.upperBand") := NA]
        }
    }

    ## export
    return(object)
}

######################################################################
### confint.ate.R ends here
