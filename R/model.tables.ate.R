### model.tables.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 16 2024 (11:48) 
## Version: 
## Last-Updated: Oct 16 2024 (12:47) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * model.tables.ate (documentation)
##' @title Statistical Inference for the Average Treatment Effect
##' @description Export estimated average treatment effects with their uncertainty (standard errors, confidence intervals and p-values).
##'
##' @param object A \code{ate} object, i.e. output of the \code{ate} function.
##' @param contrasts [character vector] levels of the treatment variable for which the estimates should be assessed or compared. Default is to consider all levels.
##' @param times [numeric vector] The timepoints at which the estimates should be displayed. Default is to consider all timepoints.
##' @param estimator [character] The type of estimator relative to which the estimates should be displayed. 
##' @param type [character] should the average risk per treatment be displayed (\code{"meanRisk"}),
##' or the difference in average risk between any two pairs of treatments (\code{"diffRisk"}),
##' or the ratio in average risk between any two pairs of treatments (\code{"ratioRisk"}).
##' @param ... Additional arguments (meanRisk.transform, ..., method.band, ...) passed to \code{\link{confint.ate}}.
##' 
##' @return a data.frame.
##' 
##' @author Brice Ozenne \email{broz@@sund.ku.dk}


## * model.tables.ate (code)
##' @export
model.tables.ate <- function(object, contrasts = NULL, times = NULL, estimator = NULL, type = NULL, ...){

    ## *** check arguments
    if(is.null(estimator)){
        estimator <- object$estimator[1]
    }else{
        estimator <- match.arg(estimator, object$estimator)
    }
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
    if(is.null(times)){
        times <- object$eval.times
    }else{
        if(any(times %in% object$eval.times == FALSE)){
            stop("Unknown timepoint ",paste(setdiff(times,object$eval.times), collapse = ", ")," in argument \'times\'. \n",
                 "Valid timepoints: ",paste(object$eval.times, collapse = ", "),". \n")
        }
    }
    if(is.null(type)){
        type <- "meanRisk"
    }else{
        type <- match.arg(type, c("meanRisk","diffRisk","ratioRisk"))
    }

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(object$inference$se == FALSE){
        stop("Cannot evaluate the uncertainty about the estimates when the standard error has not been stored. \n",
             "Set argument \'se\' to TRUE when calling the ate function \n")
    }
    
    ## *** reduce object
    object.reduce <- object
    subset.meanRisk <- which(object$meanRisk$estimator %in% estimator & object$meanRisk$time %in% times & object$meanRisk$treatment %in% contrasts)
    object.reduce$meanRisk <- object$meanRisk[subset.meanRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$meanRisk$treatment)){
        object.reduce$meanRisk$treatment <- droplevels(object.reduce$meanRisk$treatment)
    }

    subset.diffRisk <- which(object$diffRisk$estimator %in% estimator & object$diffRisk$time %in% times & object$diffRisk$B %in% contrasts & object$diffRisk$A %in% contrasts)
    object.reduce$diffRisk <- object$diffRisk[subset.diffRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$diffRisk$A)){
        object.reduce$diffRisk$A <- droplevels(object.reduce$diffRisk$A)
    }
    if(is.factor(object.reduce$diffRisk$B)){
        object.reduce$diffRisk$B <- droplevels(object.reduce$diffRisk$B)
    }

    subset.ratioRisk <- which(object$ratioRisk$estimator %in% estimator & object$ratioRisk$time %in% times & object$ratioRisk$B %in% contrasts & object$ratioRisk$A %in% contrasts)
    object.reduce$ratioRisk <- object$ratioRisk[subset.ratioRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$ratioRisk$A)){
        object.reduce$ratioRisk$A <- droplevels(object.reduce$ratioRisk$A)
    }
    if(is.factor(object.reduce$ratioRisk$B)){
        object.reduce$ratioRisk$B <- droplevels(object.reduce$ratioRisk$B)
    }

    object.reduce$iid <- list(lapply(object.reduce$iid[[estimator]][contrasts], function(iIID){iIID[,object$eval.times %in% times, drop = FALSE]}))
    names(object.reduce$iid) <- estimator
    object.reduce$estimator <- estimator ## side effect: drop attributes but they are not used by confintIID.ate
    object.reduce$eval.times <- times
    object.reduce$contrasts <- contrasts
    object.reduce$allContrasts <- utils::combn(contrasts, m = 2)
    object.reduce$inference.contrasts <- contrasts
    object.reduce$inference.allContrasts <- utils::combn(contrasts, m = 2)

    ## *** call confint
    out <- stats::confint(object.reduce)[[type]]

    ## *** export
    return(out)

}


##----------------------------------------------------------------------
### model.tables.ate.R ends here
