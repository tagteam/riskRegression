### coxCenter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:31) 
## Version: 
## Last-Updated: feb 16 2026 (09:48) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# {{{ coxCenter
## * coxCenter
#' @title Extract the mean value of the covariates
#' @description Extract the mean value of the covariates
#' @name coxCenter 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxCenter
#' @export
coxCenter <- function(object){
  UseMethod("coxCenter") 
} 

## ** coxCenter.cph
#' @rdname coxCenter
#' @method coxCenter cph
#' @export
coxCenter.cph <- function(object){
  return(setNames(object$means, object$mmcolnames))
}

## ** coxCenter.coxph
#' @rdname coxCenter
#' @method coxCenter coxph
#' @export
coxCenter.coxph <- function(object){
  return(setNames(object$means, names(coef(object))))
}


## ** coxCenter.phreg
#' @rdname coxCenter
#' @method coxCenter phreg
#' @export
coxCenter.phreg <- function(object){
    return(apply(object$X,2,mean))
}

# }}}

######################################################################
### coxCenter.R ends here
