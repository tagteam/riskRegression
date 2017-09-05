### predictCoxPL.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (16:43) 
## Version: 
## last-updated: sep  5 2017 (11:02) 
##           By: Brice Ozenne
##     Update #: 57
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Computation of survival probabilities from Cox regression models using the product limit estimator.
#' @description Same as predictCox except that the survival is estimated using the product limit estimator.
#' @examples 
#' library(survival)
#' 
#' set.seed(10)
#' d <- sampleData(40,outcome="survival")
#' nd <- sampleData(4,outcome="survival")
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#' predictCoxPL(fit, newdata = d, times = 1:5)
#' fit <- coxph(Surv(time,event)~X1 + X2 + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#' predictCoxPL(fit, newdata = d, times = 1:5)
#' 
#'
#' #### Compare exp to product limit
#' set.seed(10)
#' A <- predictCoxPL(fit, newdata = d[1:5], times = 1:5, se = TRUE, band = TRUE, logTransform = FALSE)
#' set.seed(10)
#' B <- predictCox(fit, newdata = d[1:5], times = 1:5, se = TRUE, band = TRUE, logTransform = FALSE)
#'
#' A$survival - B$survival
#' A$survival.lower - B$survival.lower
#' A$survival.upper - B$survival.upper
#' A$survival.lowerBand - B$survival.lowerBand
#' A$survival.upperBand - B$survival.upperBand
#' @export
predictCoxPL <- function(object,
                         newdata,
                         times,
                         type = c("cumhazard","survival"),
                         se = FALSE,
                         band = FALSE,
                         ...){

    # {{{ normalize arguments
    if(is.data.table(newdata)){
        newdata <- copy(newdata)
    }else{
        newdata <- as.data.table(newdata)
    }
    
    if("survival" %in% type == FALSE){
        stop("The argument \'type\' must contain \"survival\" \n",
             "use the function predictCox otherwise \n")
    }
    # }}}

    # {{{ run original prediction
    original.res <- predictCox(object = object,
                               newdata = newdata,
                               times = times,
                               type = type,
                               se = se,
                               band = band,
                               ...)
    logTransform <- class(original.res$transformation.survival)=="function"
    infoVar <- CoxVariableName(object)
    X.design <- as.data.table(CoxDesign(object))
    # }}}
    
    # {{{ compute survival
    if(infoVar$is.strata){

        object.strata <- CoxStrata(object, data = NULL, stratavars = infoVar$stratavars)
        object.levelStrata <- levels(object.strata)
        new.strata <- CoxStrata(object, data = newdata, 
                                sterms = infoVar$sterms, 
                                stratavars = infoVar$stratavars, 
                                levels = object.levelStrata, 
                                stratalevels = infoVar$stratalevels)

        Ustrata <- unique(new.strata)
        n.Ustrata <- length(Ustrata)

        for(iStrata in 1:n.Ustrata){ # iStrata <- 1
            indexStrata.object <- which(object.strata==Ustrata[iStrata])
            indexStrata.newdata <- which(new.strata==Ustrata[iStrata])
                
            all.times <- X.design[indexStrata.object,.SD$stop]
            all.times <- sort(unique(all.times[all.times <= max(times)]))
            if(length(all.times)>0){
                res.tempo <- predictCox(object,
                                        newdata = newdata[indexStrata.newdata,],
                                        times = all.times,
                                        type = "hazard")$hazard
                survival.PL <- t(apply(res.tempo, 1, function(x){cumprod(1-x)}))
                index.jump <- prodlim::sindex(eval.times = times, jump.times = all.times)
                original.res$survival[indexStrata.newdata,] <-  survival.PL[,index.jump,drop=FALSE]
            }
        }
        
    }else{
        all.times <- unique(sort(X.design$stop[X.design$stop <= max(times)]))
        if(length(all.times)>0){
            res.tempo <- predictCox(object,
                                    newdata = newdata,
                                    times = all.times,
                                    type = "hazard")$hazard
            survival.PL <- t(apply(res.tempo, 1, function(x){cumprod(1-x)}))
            index.jump <- prodlim::sindex(eval.times = times, jump.times = all.times)
            original.res$survival <-  survival.PL[,index.jump,drop=FALSE]
        }
    }
    # }}}
    
    # {{{ update confidence intervals
    if(se){
        zval <- qnorm(1- (1-original.res$conf.level)/2, 0,1)
        if(logTransform){
            original.res$survival.lower <- exp(-exp(log(-log(original.res$survival)) + zval*original.res$survival.se))
            original.res$survival.upper <- exp(-exp(log(-log(original.res$survival)) - zval*original.res$survival.se))                
        }else{
            # to keep matrix format even when out$survival contains only one line
            original.res$survival.lower <- original.res$survival.upper <- matrix(NA,
                                                                                 nrow = NROW(original.res$survival.se),
                                                                                 ncol = NCOL(original.res$survival.se)) 
            original.res$survival.lower[] <- apply(original.res$survival - zval*original.res$survival.se,2,pmax,0)
            original.res$survival.upper[] <- apply(original.res$survival + zval*original.res$survival.se,2,pmin,1)                        
        }
    }
    # }}}

    # {{{ update confidence bands
    if(band){
        quantile95 <- colMultiply_cpp(original.res$survival.se,original.res$quantile.band)

        if(logTransform){
            original.res$survival.lowerBand <- exp(-exp(log(-log(original.res$survival)) + quantile95))
            original.res$survival.upperBand <- exp(-exp(log(-log(original.res$survival)) - quantile95))
        }else{
            original.res$survival.lowerBand <- original.res$survival.upperBand <- matrix(NA, nrow = NROW(original.res$survival.se), ncol = NCOL(original.res$survival.se)) 
            original.res$survival.lowerBand[] <- apply(original.res$survival - quantile95,2,pmax,0)
            original.res$survival.upperBand[] <- apply(original.res$survival + quantile95,2,pmin,1)
        }
        
    }
    # }}}
    
    return(original.res)
    
    
}
#----------------------------------------------------------------------
### predictCoxPL.R ends here
