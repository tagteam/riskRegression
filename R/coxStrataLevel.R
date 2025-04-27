### coxStrataLevel.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:35) 
## Version: 
## Last-Updated: Apr 27 2025 (07:35) 
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

#' @title Returns the name of the strata in Cox model
#' @description Return the name of the strata in Cox model
#' @name coxStrataLevel
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#'
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxStrataLevel
#' @export
coxStrataLevel <- function(object) UseMethod("coxStrataLevel")

## ** coxStrataLevel.coxph
#' @rdname coxStrataLevel
#' @method coxStrataLevel coxph
#' @export
coxStrataLevel.coxph <- function(object){
    if(!is.null(object$strata)){
        return(levels(object$strata))
    }else{
        return(NULL)
    }
}

## ** coxStrataLevel.cph
#' @rdname coxStrataLevel
#' @method coxStrataLevel cph
#' @export
coxStrataLevel.cph <- function(object){
    if(!is.null(object$strata)){
        return(levels(object$strata))
    }else{
        return(NULL)
    }
}

## ** coxStrataLevel.phreg
#' @rdname coxStrataLevel
#' @method coxStrataLevel phreg
#' @export
coxStrataLevel.phreg <- function(object){
    if(!is.null(object$strata.name)){
       return( levels(object$model.frame[[object$strata.name]]) )
    }else{
        return( NULL )
    }
}

## ** coxStrataLevel.prodlim
#' @rdname coxStrataLevel
#' @method coxStrataLevel prodlim
#' @export
coxStrataLevel.prodlim <- function(object){
    if(!is.null(object$xlevels)){
       return( levels(interaction(expand.grid(object$xlevels[names(object$X)]),sep = ", ")) )
    }else{
        return( NULL )
    }
}


######################################################################
### coxStrataLevel.R ends here
