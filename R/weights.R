### weights.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 15 2026 (17:19) 
## Version: 
## Last-Updated: apr 17 2026 (17:36) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * weights.ate
##' @title Extract Inverse Probability Weights
##' @description Extract weights used to re-balance treatment groups (IPTW) and handle right-censoring (IPCW).
##'
##' @param object A \code{ate} object, i.e. output of the \code{ate} function.
##' @param type [character] type of weight to be extracted: \itemize{
##' \item \code{"probaT"}: probability that a given subject (in rows) get a given treatment (in columns).
##' \item \code{"probaC"}: probability that a given subject (in rows) is still at risk at a given time (in columns).
##' \item \code{"IPTW"}: observed treatment indicator divided by probaT.
##' \item \code{"IPCW"}: observed at risk indicator divided by probaT.
##' \item \code{"IPTWbox"}: create an object containing the weights and corresponding influence function (if exported by ate) that can be passed to the argument \code{treatment} of the \code{ate} function.
##' }
##' @param ... Not used. For compatibility with the generic method.
##'
##' @return A matrix with as many lines as observations and as many columns as treatment groups (\code{"probaT"} or \code{"IPTW"})
##' or as timepoint (\code{"probaC"} or \code{"IPCW"}) at which to evaluate average treatment effects.
##'
##' @details Relevant when using IPTW or AIPTW estimators but not for G-formula estimators as no weights is involved.

##' @export
weights.ate <- function(object, type = "IPTW", ...){

    ## ** normalize user input
    type <- match.arg(type,c("probaT","probaC","IPTW","IPCW","IPTWbox"))
    if(attr(object$estimator,"IPTW")==FALSE && attr(object$estimator,"AIPTW") == FALSE){
        stop("No weight to extract for Gformula estimator. \n",
             "Consider providing a treatment model, possibly also a censoring model when calling ate. \n",
             "If already done, set the argument \'estimator\' to \"IPTW\" or \"AIPTW\". \n")
    }
    if(is.null(object$weights$IPTW)){
        stop("No weight has been stored in the object. \n",
             "Consider setting the argument \'store\' to c(weights = TRUE) to store the weights. \n")
    }
    
    ## ** extract
    if(type == "IPTWbox"){
        proba <- attr(object$weights$IPTW,"proba")
        colnames(proba) <- object$contrasts

        iid.proba <- attr(object$weights$IPTW,"iid.proba")
        if(!is.null(iid.proba)){
            names(iid.proba) <- object$contrasts
        }
        out <- IPWbox(variable = object$variables["treatment"],
                      id = NULL,
                      proba = proba,
                      IF = iid.proba)

    }else if(type %in% c("IPTW","IPCW")){
        out <- object$weights[[type]]
        attr(out,"proba") <- NULL
        attr(out,"iid.proba") <- NULL
    }else if(type == "probaT"){
        out <- attr(object$weights$IPTW, "proba")
    }else if(type == "probaC"){
        out <- attr(object$weights$IPCW, "proba")
    }

    ## ** export
    return(out)
    
}

##----------------------------------------------------------------------
### weights.R ends here
