### coxVarCov.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:36) 
## Version: 
## Last-Updated: Apr 27 2025 (07:36) 
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
#' @title Extract the variance covariance matrix of the beta from a Cox model
#' @description Extract the variance covariance matrix of the beta from a Cox model
#' @name coxVarCov 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details Should return \code{NULL} if the Cox model has no covariate. 
#' The rows and columns of the variance covariance matrix must be named with the names used in the design matrix.

#' @rdname coxVarCov
#' @export
coxVarCov <- function(object){
  UseMethod("coxVarCov") 
} 

## ** coxVarCov.cph
#' @rdname coxVarCov
#' @method coxVarCov cph
#' @export
coxVarCov.cph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    colnames(Sigma) <- object$mmcolnames
    rownames(Sigma) <- object$mmcolnames
  }
  
  return(Sigma)
}

## ** coxVarCov.coxph
#' @rdname coxVarCov
#' @method coxVarCov coxph
#' @export
coxVarCov.coxph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}

## ** coxVarCov.phreg
#' @rdname coxVarCov
#' @method coxVarCov phreg
#' @export
coxVarCov.phreg <- function(object){
  
  Sigma <- -solve(object$hessian)
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}
# }}}


## ** coxVarCov.prodlim
#' @rdname coxVarCov
#' @method coxVarCov prodlim
#' @export
coxVarCov.prodlim <- function(object){
    return(NULL)
}
# }}}



######################################################################
### coxVarCov.R ends here
