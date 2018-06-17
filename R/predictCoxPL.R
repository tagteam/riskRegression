### predictCoxPL.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (16:43) 
## Version: 
## last-updated: jun 17 2018 (11:17) 
##           By: Brice Ozenne
##     Update #: 92
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * predictCoxPL (documentation)
#' @title Computation of survival probabilities from Cox regression models using the product limit estimator.
#' @description Same as predictCox except that the survival is estimated using the product limit estimator.
#' @name predictCoxPL
#'
#' @inheritParams predictCox
#' @param ... additional arguments to be passed to \code{\link{predictCox}}.
#' 
#' @examples 
#' library(survival)
#'
#' #### generate data ####
#' set.seed(10)
#' d <- sampleData(40,outcome="survival")
#' nd <- sampleData(4,outcome="survival")
#' d$time <- round(d$time,1)
#'
#' #### Cox model ####
#' fit <- coxph(Surv(time,event)~ X1 + X2 + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#'
#' ## exponential approximation
#' predictCox(fit, newdata = d, times = 1:5)
#' 
#' ## product limit
#' predictCoxPL(fit, newdata = d, times = 1:5)
#'
#' #### stratified Cox model ####
#' fitS <- coxph(Surv(time,event)~ X1 + strata(X2) + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#'
#' ## exponential approximation
#' predictCox(fitS, newdata = d, times = 1:5)
#' 
#' ## product limit
#' predictCoxPL(fitS, newdata = d, times = 1:5)
#'

## * predictCoxPL (code)
#' @rdname predictCoxPL 
#' @export
predictCoxPL <- function(object,
                         newdata,
                         times,
                         type = c("cumhazard","survival"),
                         keep.strata = TRUE,
                         keep.infoVar = FALSE,
                         ...){
    stop <- NULL ## [:CRANcheck:] data.table
    
    ## ** normalize arguments
    object.modelFrame <- coxModelFrame(object)
    if(is.data.table(newdata)){
        newdata <- copy(newdata)
    }else{
        newdata <- as.data.table(newdata)
    }
    
    if("survival" %in% type == FALSE){
        stop("The argument \'type\' must contain \"survival\" \n",
             "use the function predictCox otherwise \n")
    }

    ## ** run original prediction
    original.res <- predictCox(object = object,
                               newdata = newdata,
                               times = times,
                               type = type,
                               keep.strata = TRUE,
                               keep.infoVar = TRUE,
                               ...)
    infoVar <- original.res$infoVar
    if(keep.infoVar==FALSE){
        original.res$infoVar <- NULL
    }
    
    ## ** compute survival
    if(infoVar$is.strata){

        object.strata <- coxStrata(object, data = NULL, strata.vars = infoVar$stratavars)
        new.strata <- original.res$strata
        if(keep.strata == FALSE){
            original.res$strata <- NULL
        }
        Ustrata <- unique(new.strata)
        n.Ustrata <- length(Ustrata)

        for(iStrata in 1:n.Ustrata){ # iStrata <- 1
            indexStrata.object <- which(object.strata==Ustrata[iStrata])
            indexStrata.newdata <- which(new.strata==Ustrata[iStrata])
                
            all.times <- object.modelFrame[indexStrata.object,.SD$stop]
            all.times <- sort(unique(all.times[all.times <= max(times)]))
            if(length(all.times)>0){
                res.tempo <- predictCox(object,
                                        newdata = newdata[indexStrata.newdata,],
                                        times = all.times,
                                        type = "hazard")$hazard
                survival.PL <- cbind(1,t(apply(res.tempo, 1, function(x){cumprod(1-x)})))
                index.jump <- prodlim::sindex(eval.times = times, jump.times = all.times)+1
                original.res$survival[indexStrata.newdata,] <-  survival.PL[,index.jump,drop=FALSE]
            }
        }
        
    }else{
        all.times <- unique(sort(object.modelFrame[stop <= max(times), stop]))
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

    ## ** export
    return(original.res)
    
    
}
#----------------------------------------------------------------------
### predictCoxPL.R ends here
