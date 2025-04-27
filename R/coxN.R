### coxN.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:34) 
## Version: 
## Last-Updated: Apr 27 2025 (07:34) 
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
## * coxN
#' @title Extract the number of observations from a Cox model
#' @description Extract the number of observations from a Cox model
#' @name coxN 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#'

#' @rdname coxN
#' @export
coxN <- function(object){
  UseMethod("coxN") 
} 

## ** coxN.cph
#' @rdname coxN
#' @method coxN cph
#' @export
coxN.cph <- function(object){
    return(sum(object$n))
}

## ** coxN.coxph
#' @rdname coxN
#' @method coxN coxph
#' @export
coxN.coxph <- function(object){
  return(object$n)
}

## ** coxN.phreg
#' @rdname coxN
#' @method coxN default
#' @export
coxN.default <- function(object){
  return(stats::nobs(object))
}

#' @rdname coxN
#' @method coxN phreg
#' @export
coxN.phreg <- function(object){
  return(NROW(object$model.frame))
}

## ** coxN.CSC
#' @rdname coxN
#' @method coxN CauseSpecificCox
#' @export
coxN.CauseSpecificCox <- function(object){
  return(sapply(object$models,coxN))
}

## ** coxN.CSC
#' @rdname coxN
#' @method coxN glm
#' @export
coxN.glm <- function(object){
  return(stats::nobs(object))
}

## ** coxN.prodlim
#' @rdname coxN
#' @method coxN prodlim
#' @export
coxN.prodlim <- function(object){
  return(NROW(object$model.response))
}



######################################################################
### coxN.R ends here
