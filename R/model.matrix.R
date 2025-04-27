### model.matrix.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:39) 
## Version: 
## Last-Updated: Apr 27 2025 (07:40) 
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
#' @title Extract design matrix for cph objects
#' @description Extract design matrix for cph objects
#' @param object a cph object.
#' @param data a dataset.
#' @param ... not used
#' 
#' @method model.matrix cph
model.matrix.cph <- function(object, data, ...){

    if(all(object$Design$assume %in% c("category","asis","interaction","strata"))){
        out <- predict(object, newdata = data, type = "x")

        ## put back rms names
        colnames(out) <- object$Design$colnames[match(object$Design$mmcolnames,colnames(out))]
    }else{
        type.special <- setdiff(object$Design$assume,c("category","asis","interaction","strata"))
        name.special <- attr(object$terms,"term.labels")[object$Design$assume %in% type.special]
        ## warning("Special operator(s) in the formula: ",paste(name.special,collapse = ", "),".\n",
        ##         "Associated type: \"",paste(type.special,collapse = "\" \""),"\". \n",
        ##         "If it is a spline term, consider using survival::coxph instead of rms::cph. \n")

        ## set categorical variables to their appropriate level
        if(any(object$Design$assume=="category")){ 
            for(iVar in object$Design$name[object$Design$assume=="category"]){
                if(!is.factor(data[[iVar]]) || any(sort(levels(data[[iVar]])) != sort(object$Design$parms[[iVar]])) ){
                    data[[iVar]] <- factor(data[[iVar]], levels=object$Design$parms[[iVar]])
                }
            }
        }

        ## borrowed from survival::coxph
        mf <- stats::model.frame(object)
        Terms <- stats::delete.response(terms(mf))
        out <- model.matrix(Terms, data)

        ## remove intercept
        out <- out[,attr(out,"assign")!=0,drop=FALSE]

        ## put back rms names
        if(!is.null(attr(object$Design$mmcolnames,"alt"))){
            colnames(out) <- object$Design$colnames[match(attr(object$Design$mmcolnames,"alt"),colnames(out))]
        }else{
            colnames(out) <- object$Design$colnames[match(object$Design$mmcolnames,colnames(out))]
        }
        ## add strata if possible
        try(attr(out,"strata") <- attr(predict(object, newdata = data, type = "x"),"strata"), silent = TRUE)
    }
    return(out)
    
}

## ** model.matrix.phreg
#' @title Extract design matrix for phreg objects
#' @description Extract design matrix for phreg objects
#' @param object a phreg object.
#' @param data a dataset.
#' @param ... not used
#' 
#' @details mainly a copy paste of the begining of the \code{phreg} function.
#' 
#' @method model.matrix phreg
model.matrix.phreg <- function(object, data, ...){
    special <- c("strata", "cluster")
    Terms <- stats::delete.response(attr(object$model.frame, "terms"))

    ## remove specials
    if (length(attributes(Terms)$specials$cluster)>0) {
        ts <- survival::untangle.specials(Terms, "cluster")
        Terms <- Terms[-ts$terms]
    }
    if (length(attributes(Terms)$specials$strata)>0) {
        ts <- survival::untangle.specials(Terms, "strata")
        Terms <- Terms[-ts$terms]
    }
    attr(Terms,"intercept") <- 1 ## keep intercept to have the same behavior with and without categorical variables 
    X <- model.matrix(Terms, data)
    return(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE])
}



######################################################################
### model.matrix.R ends here
