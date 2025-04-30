### coxFormula.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:31) 
## Version: 
## Last-Updated: Apr 29 2025 (06:51) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * coxFormula
#' @title Extract the formula from a Cox model
#' @description Extract the formula from a Cox model
#' @name coxFormula 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxFormula
#' @export
coxFormula <- function(object){
  UseMethod("coxFormula") 
} 

## ** coxFormula.cph
#' @rdname coxFormula
#' @method coxFormula cph
#' @export
coxFormula.cph <- function(object){
  return(object$sformula)
}

## ** coxFormula.coxph
#' @rdname coxFormula
#' @method coxFormula coxph
#' @export
coxFormula.coxph <- function(object){
    if(object$nevent>0){
        return(object$formula)
    }else{
        out <- object$terms
        extra.attr <- names(attributes(out))
        for(iAttr in extra.attr){
            attr(out, iAttr) <- NULL
        }
        return(as.formula(out))
    }
}

## ** coxFormula.phreg
#' @rdname coxFormula
#' @method coxFormula phreg
#' @export
coxFormula.phreg <- function(object){
  return(object$formula)
}
## ** coxFormula.glm
#' @rdname coxFormula
#' @method coxFormula glm
#' @export
coxFormula.glm <- function(object){
  return(stats::formula(object))
}

## ** coxFormula.prodlim
#' @rdname coxFormula
#' @method coxFormula prodlim
#' @export
coxFormula.prodlim <- function(object){
    return(stats::formula(object))
}

## ** coxFormula.coxnet
#' @rdname coxFormula
#' @method coxFormula coxnet
#' @export
coxFormula.coxnet <- function(object){
    return(object$formula)
}

######################################################################
### coxFormula.R ends here
