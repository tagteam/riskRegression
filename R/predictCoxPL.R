### predictCoxPL.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (16:43) 
## Version: 
## last-updated: maj 28 2018 (17:56) 
##           By: Brice Ozenne
##     Update #: 66
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
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata A \code{data.frame} or \code{data.table} containing
#'     the values of the predictor variables defining subject specific
#'     predictions. Should have the same structure as the data set
#'     used to fit the \code{object}.
#' @param times Time points at which to evaluate the predictions.
#' @param type the type of predicted value. 
#'     Choices are \code{"hazard"}, \code{"cumhazard"}, and \code{"survival"}. 
#'     See \code{\link{predictCox}} for more details.
#' @param se Logical. If \code{TRUE} add the standard error to the output.
#' @param band Logical. If \code{TRUE} add the confidence band to the output.
#' @param ... additional arguments to be passed to \code{\link{predictCox}}.
#' 
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
#' A <- predictCoxPL(fit, newdata = d[1:5], times = 1:5, se = TRUE, band = TRUE)
#' set.seed(10)
#' B <- predictCox(fit, newdata = d[1:5], times = 1:5, se = TRUE, band = TRUE)
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
    infoVar <- coxVariableName(object)
    X.design <- as.data.table(coxDesign(object))
    # }}}
    
    # {{{ compute survival
    if(infoVar$is.strata){

        object.strata <- coxStrata(object, data = NULL, strata.vars = infoVar$strata.vars)
        object.levelStrata <- levels(object.strata)
        new.strata <- coxStrata(object, data = newdata, 
                                sterms = infoVar$sterms, 
                                strata.vars = infoVar$strata.vars, 
                                levels = object.levelStrata, 
                                strata.levels = infoVar$strata.levels)

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
    
    return(original.res)
    
    
}
#----------------------------------------------------------------------
### predictCoxPL.R ends here
