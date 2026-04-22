### model.tables.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 16 2024 (11:48) 
## Version: 
## Last-Updated: apr 22 2026 (18:03) 
##           By: Brice Ozenne
##     Update #: 37
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
##' @param x A \code{ate} object, i.e. output of the \code{ate} function.
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
model.tables.ate <- function(x, contrasts = NULL, times = NULL, estimator = NULL, type = NULL, ...){

    ## **  normalize user input
    ## *** estimator
    ## rename estimator --> user.estimator to avoid confusion with the column names
    if(is.null(estimator)){
        user.estimator <- x$estimator[1]
    }else{
        estimator <- toupper(estimator)
        if(length(estimator) == 1 && estimator == "ALL"){
            user.estimator <- x$estimator
        }else{
            user.estimator <- match.arg(estimator, x$estimator, several.ok = TRUE)
        }
    }

    ## *** contrasts
    if(is.null(contrasts)){
        contrasts <- x$contrasts
    }else{
        contrasts.original <- contrasts
        contrasts <- match.arg(contrasts, x$contrasts, several.ok = TRUE)
        if(length(contrasts.original)!=length(contrasts)){
            stop("Unknown value ",paste(setdiff(contrasts.original,x$contrasts), collapse = ", ")," in argument \'contrasts\'. \n",
                 "Valid values: ",paste(x$contrasts, collapse = ", "),". \n")
        }
    }

    ## *** times
    if(is.null(times)){
        times <- x$eval.times
    }else{
        if(any(times %in% x$eval.times == FALSE)){
            stop("Unknown timepoint ",paste(setdiff(times,x$eval.times), collapse = ", ")," in argument \'times\'. \n",
                 "Valid timepoints: ",paste(x$eval.times, collapse = ", "),". \n")
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

    ## ** reduce object
    object.reduce <- x
    subset.meanRisk <- which(x$meanRisk$estimator %in% user.estimator & x$meanRisk$time %in% times & x$meanRisk$treatment %in% contrasts)
    object.reduce$meanRisk <- x$meanRisk[subset.meanRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$meanRisk$treatment)){
        object.reduce$meanRisk$treatment <- droplevels(object.reduce$meanRisk$treatment)
    }

    subset.diffRisk <- which(x$diffRisk$estimator %in% user.estimator & x$diffRisk$time %in% times & x$diffRisk$B %in% contrasts & x$diffRisk$A %in% contrasts)
    object.reduce$diffRisk <- x$diffRisk[subset.diffRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$diffRisk$A)){
        object.reduce$diffRisk$A <- droplevels(object.reduce$diffRisk$A)
    }
    if(is.factor(object.reduce$diffRisk$B)){
        object.reduce$diffRisk$B <- droplevels(object.reduce$diffRisk$B)
    }

    subset.ratioRisk <- which(x$ratioRisk$estimator %in% user.estimator & x$ratioRisk$time %in% times & x$ratioRisk$B %in% contrasts & x$ratioRisk$A %in% contrasts)
    object.reduce$ratioRisk <- x$ratioRisk[subset.ratioRisk] ## does not work directly due to confusion from data.table between values and column names (e.g. estimator is both)
    if(is.factor(object.reduce$ratioRisk$A)){
        object.reduce$ratioRisk$A <- droplevels(object.reduce$ratioRisk$A)
    }
    if(is.factor(object.reduce$ratioRisk$B)){
        object.reduce$ratioRisk$B <- droplevels(object.reduce$ratioRisk$B)
    }
    if(x$inference$iid){
        object.reduce$iid <- lapply(user.estimator, function(iE){stats::setNames(lapply(contrasts, function(iC){object.reduce$iid[[iE]][[iC]][,x$eval.times %in% times, drop = FALSE]}), contrasts)})
        names(object.reduce$iid) <- user.estimator
    }else if(x$inference$bootstrap){
        object.reduce$boot$t0 <- x$boot$t0[c(subset.meanRisk,NROW(x$meanRisk)+subset.diffRisk,NROW(x$meanRisk)+NROW(x$diffRisk)+subset.ratioRisk)]
        object.reduce$boot$t <- x$boot$t[,c(subset.meanRisk,NROW(x$meanRisk)+subset.diffRisk,NROW(x$meanRisk)+NROW(x$diffRisk)+subset.ratioRisk),drop=FALSE]
    }
    object.reduce$estimator <- user.estimator ## side effect: drop attributes but they are not used by confintIID.ate
    object.reduce$eval.times <- times
    object.reduce$contrasts <- contrasts
    object.reduce$allContrasts <- utils::combn(contrasts, m = 2)
    object.reduce$inference.contrasts <- contrasts
    object.reduce$inference.allContrasts <- utils::combn(contrasts, m = 2)

    ## *** call confint
    out <- stats::confint(object.reduce, ...)[[type]]

    ## *** export
    return(out[])

}


##----------------------------------------------------------------------
### model.tables.ate.R ends here
