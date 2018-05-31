### confint.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:08) 
## Version: 
## Last-Updated: maj 31 2018 (11:49) 
##           By: Brice Ozenne
##     Update #: 298
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
##' @param conf.level Level of confidence.
##' @param transform.absRisk the transformation used to improve coverage
##' of the confidence intervals for the predicted absolute risk in small samples.
##' @param nsim.band the number of simulations used to compute the quantiles for the confidence bands.
##' @param seed Integer passed to set.seed when performing simulation for the confidence bands.
##' If not given or NA no seed is set. 
##' @param ... not used.
##'
##' @details The confidence bands and confidence intervals are automatically restricted to the interval [0;1].
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

## * confint.ate (code)
##' @rdname confint.ate
##' @export
confint.ate <- function(object,
                        conf.level = 0.95,
                        nsim.band = 1e4,
                        ate.transform = "none",
                        diffAte.transform = "none",
                        ratioAte.transform = "none",
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
                                      conf.level = conf.level)

        }else{
            object <- confintIID.ate(object,
                                     nsim.band = nsim.band,
                                     ate.transform = ate.transform,
                                     diffAte.transform = diffAte.transform,
                                     ratioAte.transform = ratioAte.transform,
                                     seed = seed,
                                     conf.level = conf.level)
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
                           ate.transform,
                           diffAte.transform,
                           ratioAte.transform,
                           seed){
    
    ## ** check arguments
    if(!is.null(object$ate.transform) && object$ate.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    if(!is.null(object$diffAte.transform) && object$diffAte.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    if(!is.null(object$diffRatio.transform) && object$diffRatio.transform != "none"){
        stop("Cannot work with standard errors that have already been transformed \n")
    }
    object$ate.transform <- match.arg(ate.transform, c("none","log","loglog","cloglog"))
    object$diffAte.transform <- match.arg(diffAte.transform, c("none","atanh"))
    object$ratioAte.transform <- match.arg(ratioAte.transform, c("none","log"))

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
    if(object$ate.transform != "none"){
        ## transform standard error
        object$meanRisk[,c("se") := transformSE(estimate = object$meanRisk[["meanRisk"]],
                                                se = object$ate.se,
                                                type = object$ate.transform)[,1]]

        if(object$band){
            object$ate.iid <- transformIID(estimate = object$meanRisk[["meanRisk"]],
                                           iid = object$ate.iid,
                                           type = object$ate.transform,
                                           format = "matrix")
        }

    }else{
        object$meanRisk[,c("se") := object$ate.se]
    } 
    object$ate.se <- NULL ## remove untransformed se to avoid confusion
    
    ate.min.value <- switch(object$ate.transform,
                            "none" = 0,
                            "log" = NULL,
                            "loglog" = NULL,
                            "cloglog" = NULL)
    ate.max.value <- switch(object$ate.transform,
                            "none" = 1,
                            "log" = 1,
                            "loglog" = NULL,
                            "cloglog" = NULL)
    ## confidence intervals
    if(object$se){
        object$meanRisk[,c("lower","upper") := transformCI(estimate = object$meanRisk[["meanRisk"]],
                                                           se = object$meanRisk[["se"]],
                                                           quantile = zval,
                                                           type = object$ate.transform,
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
                                                                   type = object$ate.transform,
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
    if(object$diffAte.transform != "none"){
        ## transform standard error
        object$riskComparison[,c("diff.se") := transformSE(estimate = object$riskComparison[["diff"]],
                                                           se = object$diffAte.se,
                                                           type = object$diffAte.transform)[,1]]

    }else{
        object$riskComparison[,c("diff.se") := object$diffAte.se]
    }
    object$diffAte.se <- NULL ## remove untransformed se to avoid confusion
    diffAte.min.value <- switch(object$diffAte.transform,
                                "none" = -1,
                                "atanh" = NULL)
    diffAte.max.value <- switch(object$diffAte.transform,
                                "none" = 1,
                                "atanh" = NULL)
    ## confidence intervals
    if(object$se){
        object$riskComparison[,c("diff.lower","diff.upper") := transformCI(estimate = object$riskComparison[["diff"]],
                                                                           se = object$riskComparison[["diff.se"]],
                                                                           quantile = zval,
                                                                           type = object$diffAte.transform,
                                                                           format = "vector",
                                                                           min.value = diffAte.min.value,
                                                                           max.value = diffAte.max.value)]
        ## p.value
        object$riskComparison[,c("diff.p.value") := transformP(estimate = object$riskComparison[["diff"]],
                                                               se = object$riskComparison[["diff.se"]],
                                                               type = object$diffAte.transform)]

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
                                                                                   type = object$diffAte.transform,
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
    if(object$ratioAte.transform != "none"){
        ## transform standard error
        object$riskComparison[,c("ratio.se") := transformSE(estimate = object$riskComparison[["ratio"]],
                                                           se = object$ratioAte.se,
                                                           type = object$ratioAte.transform)[,1]]

    }else{
        object$riskComparison[,c("ratio.se") := object$ratioAte.se]
    }
    object$ratioAte.se <- NULL ## remove untransformed se to avoid confusion
    ratioAte.min.value <- switch(object$ratioAte.transform,
                                 "none" = 0,
                                 "log" = NULL)
    ratioAte.max.value <- NULL
    
    ## confidence intervals
    if(object$se){
        object$riskComparison[,c("ratio.lower","ratio.upper") := transformCI(estimate = object$riskComparison[["ratio"]],
                                                                             se = object$riskComparison[["ratio.se"]],
                                                                             quantile = zval,
                                                                             type = object$ratioAte.transform,
                                                                             format = "vector",
                                                                             min.value = ratioAte.min.value,
                                                                             max.value = ratioAte.max.value)]
        ## p.value
        object$riskComparison[,c("ratio.p.value") := transformP(estimate = object$riskComparison[["ratio"]],
                                                                se = object$riskComparison[["ratio.se"]],
                                                                type = object$ratioAte.transform)]

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
                                                                                     type = object$ratioAte.transform,
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
