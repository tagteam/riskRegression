### vcov.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct 16 2024 (11:47) 
## Version: 
## Last-Updated: Oct 16 2024 (11:48) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov.ate (documentation)
##' @title Variance-Covariance Matrix for the Average Treatment Effect.
##' @description Variance covariance matrix for the estimated average treatment effect.
##'
##' @param object A \code{ate} object, i.e. output of the \code{ate} function.
##' @param contrasts [character vector] levels of the treatment variable for which the variance-covariance matrix should be assessed. Default is to consider all levels.
##' @param times [numeric vector] The timepoints at which the variance-covariance matrix should be displayed. Default is to consider all timepoints.
##' @param estimator [character] The type of estimator relative to which the variance-covariance matrix should be displayed. 
##' @param type [character] should the variance-covariance matrix w.r.t. the average risk per treatment be displayed (\code{"meanRisk"}),
##' or the difference in average risk between any two pairs of treatments (\code{"diffRisk"}),
##' or the ratio in average risk between any two pairs of treatments (\code{"ratioRisk"}).
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @return A numeric matrix.
##' 
##' @author Brice Ozenne \email{broz@@sund.ku.dk}

## * vcov.ate (code)
##' @export
vcov.ate <- function(object, contrasts = NULL, times = NULL, estimator = NULL, type = NULL, ...){

    ## *** normalize user input
    if(is.null(estimator)){
        estimator <- object$estimator[1]
    }else{
        estimator <- match.arg(estimator, object$estimator)
    }

    if(is.null(type)){
        type <- "meanRisk"
    }else{
        type <- match.arg(type, c("meanRisk","diffRisk","ratioRisk"))
    }

    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }

    if(is.null(object$iid)){
        stop("Missing iid decomposition in object \n",
             "Consider setting the argument \'iid\' to TRUE when calling the ate function \n")
    }

    if(object$inference$iid == FALSE){
        stop("Cannot evaluate the variance-covariance matrix without the influence function. \n",
             "Set argument \'se\' to TRUE when calling the ate function \n")
    }
    
    ## *** extract
    object.iid <- do.call(cbind,object$iid[[estimator]])
    if(type == "meanRisk"){
        out <- crossprod(object.iid)
    }else{
        n.times <- length(object$eval.times)
        vec.contrast1 <- as.numeric(factor(object$allContrasts[1,], object$contrast))
        vec.contrast2 <- as.numeric(factor(object$allContrasts[2,], object$contrast))
        
        vecTime.contrast1 <- unlist(lapply(vec.contrast1, function(iC){n.times*(iC-1)+1:n.times}))
        vecTime.contrast2 <- unlist(lapply(vec.contrast2, function(iC){n.times*(iC-1)+1:n.times}))
        
        if(type == "diffRisk"){

            out <- crossprod(object.iid[ , vecTime.contrast2, drop=FALSE] - object.iid[, vecTime.contrast1, drop=FALSE])

        }else if(type == "ratioRisk"){

            object.mean <- coef(object, estimator = estimator, type = "meanRisk")
            iid.num <- rowScale_cpp(object.iid[ , vecTime.contrast2, drop=FALSE], object.mean[vecTime.contrast1])
            iid.denum <- rowMultiply_cpp(object.iid[, vecTime.contrast1, drop=FALSE], object.mean[vecTime.contrast2]/object.mean[vecTime.contrast1]^2)
            out <- crossprod(iid.num - iid.denum)

        }
    }

    ## *** export
    name.out <- names(coef(object, estimator = estimator, type = type))
    dimnames(out) <- list(name.out,name.out)

    if(!is.null(times) || !is.null(contrasts)){
        name.subset <- names(coef(object, estimator = estimator, type = type, contrasts = contrasts, times = times))
        out <- out[name.subset,name.subset,drop=FALSE]
    }
    return(out)
}


##----------------------------------------------------------------------
### vcov.ate.R ends here
