### coxSpecialStrata.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:34) 
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

## * coxSpecial
#' @title Special characters in Cox model
#' @description Return the special character(s) of the Cox model, e.g. used to indicate the strata variables.
#' @name coxSpecial
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#'
#' @details Must return a list with at least one element strata
#' indicating the character in the formula marking the variable(s) defining the strata.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxSpecial
#' @export
coxSpecial <- function(object) UseMethod("coxSpecial")

## ** coxSpecial.coxph
#' @rdname coxSpecial
#' @method coxSpecial coxph
#' @export
coxSpecial.coxph <- function(object){
    return(list(strata = "strata",
                cluster = "cluster"))
}

## ** coxSpecial.cph
#' @rdname coxSpecial
#' @method coxSpecial cph
#' @export
coxSpecial.cph <- function(object){
  return(list(strata = "strat"))
}

## ** coxSpecial.phreg
#' @rdname coxSpecial
#' @method coxSpecial phreg
#' @export
coxSpecial.phreg <- function(object){
    return(list(strata = "strata",
                cluster = "cluster"))
}
# }}}
## ** coxSpecial.prodlim
#' @rdname coxSpecial
#' @method coxSpecial prodlim
#' @export
coxSpecial.prodlim <- function(object){
    return(NULL)
}

######################################################################
### coxSpecialStrata.R ends here
