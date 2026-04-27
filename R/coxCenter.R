### coxCenter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:31) 
## Version: 
## Last-Updated: Apr 27 2026 (12:16) 
##           By: Brice Ozenne
##     Update #: 20
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
#' @param as.factor [logical] should the original factor level be returned instead of the corresponding dummy variables? Will return a data.frame instead of a matrix as side effect.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxCenter
#' @export
coxCenter <- function(object, as.factor){
  UseMethod("coxCenter") 
} 

## ** coxCenter.cph
#' @rdname coxCenter
#' @method coxCenter cph
#' @export
coxCenter.cph <- function(object, as.factor = FALSE){
    out <- setNames(object$means, names(coef(object)))
    if(as.factor){
        out <- as.data.frame(as.list(out))
    }
  return(out)
}

## ** coxCenter.coxph
#' @rdname coxCenter
#' @method coxCenter coxph
#' @export
coxCenter.coxph <- function(object, as.factor = FALSE){
    if(as.factor & length(object$xlevels)>0){
        object.info <- coxVariableName(object, model.frame = coxModelFrame(object))

        var.factor <-  intersect(names(object$xlevels),
                                 setdiff(object.info$lpvars.original,object.info$lpvars)) ## could contain splines names e.g. ns(X,1)
        ref.factor <- sapply(object$xlevels[var.factor],"[",1)
        ref.num <- object$means[setdiff(names(object$means),object.info$lpvars[which(object.info$lpvars.original %in% var.factor)])]        
        out <- as.data.frame(c(as.list(ref.num),as.list(ref.factor)))
    }else{
        out <- setNames(object$means, names(coef(object)))
    }
    return(out)
}


## ** coxCenter.phreg
#' @rdname coxCenter
#' @method coxCenter phreg
#' @export
coxCenter.phreg <- function(object, as.factor = FALSE){

    out <- apply(object$X,2,mean)
    if(as.factor){
        out <- as.data.frame(as.list(out))
    }
    return(out)
}

# }}}

######################################################################
### coxCenter.R ends here
