### coxBaseEstimator.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:30) 
## Version: 
## Last-Updated: Apr 27 2025 (07:30) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# {{{ coxBaseEstimator
## * coxBaseEstimator
#' @title Extract the type of estimator for the baseline hazard
#' @description Extract the type of estimator for the baseline hazard
#' @name coxBaseEstimator 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxBaseEstimator
#' @export
coxBaseEstimator <- function(object){
  UseMethod("coxBaseEstimator") 
} 

## ** coxBaseEstimator.coxph
#' @rdname coxBaseEstimator
#' @method coxBaseEstimator coxph
#' @export
coxBaseEstimator.coxph <- function(object){
  return(object$method)
}

## ** coxBaseEstimator.phreg
#' @rdname coxBaseEstimator
#' @method coxBaseEstimator phreg
#' @export
coxBaseEstimator.phreg <- function(object){
  return("breslow")
}

## ** coxBaseEstimator.prodlim
#' @rdname coxBaseEstimator
#' @method coxBaseEstimator prodlim
#' @export
coxBaseEstimator.prodlim <- function(object){
    if("method.ties" %in% names(object)){
        return(object$method.ties)
    }else{
        return("breslow")
    }
}

# }}}



######################################################################
### coxBaseEstimator.R ends here
