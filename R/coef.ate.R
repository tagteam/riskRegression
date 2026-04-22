### coef.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 16 2024 (11:47) 
## Version: 
## Last-Updated: apr 22 2026 (11:24) 
##           By: Brice Ozenne
##     Update #: 17
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * coef.ate (documentation)
##' @title Estimated Average Treatment Effect.
##' @description Estimated average treatment effect.
##'
##' @param object A \code{ate} object, i.e. output of the \code{ate} function.
##' @param contrasts [character vector] levels of the treatment variable for which the estimates should be assessed or compared. Default is to consider all levels.
##' @param times [numeric vector] The timepoints at which the estimates should be displayed. Default is to consider all timepoints.
##' @param estimator [character] The type of estimator relative to which the estimates should be displayed. Can be \code{"all"} to output all estimators.
##' @param type [character] should the average risk per treatment be displayed (\code{"meanRisk"}),
##' or the difference in average risk between any two pairs of treatments (\code{"diffRisk"}),
##' or the ratio in average risk between any two pairs of treatments (\code{"ratioRisk"}).
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A numeric vector.
##' 
##' @author Brice Ozenne \email{broz@@sund.ku.dk}

## * coef.ate (code)
##' @export
coef.ate <- function(object, contrasts = NULL, times = NULL, estimator = NULL, type = NULL, ...){

    ## **  normalize user input
    ## *** estimator
    ## rename estimator --> user.estimator to avoid confusion with the column names
    if(is.null(estimator)){
        user.estimator <- object$estimator[1]
    }else{
        estimator <- toupper(estimator)
        if(length(estimator) == 1 && estimator == "ALL"){
            user.estimator <- object$estimator
        }else{
            user.estimator <- match.arg(estimator, object$estimator, several.ok = TRUE)
        }
    }


    ## *** contrasts
    if(is.null(contrasts)){
        contrasts <- object$contrasts
    }else{
        contrasts.original <- contrasts
        contrasts <- match.arg(contrasts, object$contrasts, several.ok = TRUE)
        if(length(contrasts.original)!=length(contrasts)){
            stop("Unknown value ",paste(setdiff(contrasts.original,object$contrasts), collapse = ", ")," in argument \'contrasts\'. \n",
                 "Valid values: ",paste(object$contrasts, collapse = ", "),". \n")
        }
    }

    ## *** times
    if(is.null(times)){
        times <- object$eval.times
    }else{
        if(any(times %in% object$eval.times == FALSE)){
            stop("Unknown timepoint ",paste(setdiff(times,object$eval.times), collapse = ", ")," in argument \'times\'. \n",
                 "Valid timepoints: ",paste(object$eval.times, collapse = ", "),". \n")
        }
    }

    ## *** type
    if(is.null(type)){
        type <- "meanRisk"
    }else{
        type <- match.arg(type, c("meanRisk","diffRisk","ratioRisk"))
    }

    ## *** dots
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    ## ** extract
    object.table1 <- object[[type]][object[[type]]$estimator %in% user.estimator,,drop=FALSE]
    out.times <- object.table1$time %in% times
    
    if(length(user.estimator)>1){
        name.estimator <- paste0("(t=",object.table1$time,",",object.table1$estimator,")")
    }else{
        name.estimator <- paste0("(t=",object.table1$time,")")
    }
    if(type == "meanRisk"){
        out <- setNames(object.table1$estimate, paste0(object.table1$treatment,name.estimator))
        out.contrasts <- object.table1$treatment %in% contrasts
    }else if(type == "diffRisk"){
        out <-setNames(object.table1$estimate, paste0(object.table1$B,"-",object.table1$A,name.estimator))
        out.contrasts <- (object.table1$B %in% contrasts) & (object.table1$A %in% contrasts)
    }else if(type == "ratioRisk"){
        out <- setNames(object.table1$estimate, paste0(object.table1$B,"/",object.table1$A,name.estimator))
        out.contrasts <- (object.table1$B %in% contrasts) & (object.table1$A %in% contrasts)
    }

    ## *** export
    return(out[out.times & out.contrasts])
}


##----------------------------------------------------------------------
### coef.ate.R ends here
