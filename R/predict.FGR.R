#' Predict subject specific risks (cumulative incidence) based on Fine-Gray regression model
#'
#' Predict subject specific risks (cumulative incidence) based on Fine-Gray regression model
#' @param object Result of call to \code{FGR}
#' @param newdata Predictor values of subjects for who to predict risks
#' @param times Time points at which to evaluate the risks 
#' @param ... passed to predict.crr
#'
#' @method predict FGR
#' @export
predict.FGR <- function(object,newdata,times,...){
    # {{{ check the data and the design
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (NROW(newdata) == 0) stop("No (non-missing) observations")
    rhs <- as.formula(delete.response(object$terms))
    EHF <- prodlim::model.design(rhs,
                                 newdata,
                                 dropIntercept=TRUE,
                                 specialsDesign=TRUE)
    # }}}
    # {{{ covariate design matrices
    cov1 <- cbind(EHF$cov1,EHF$design)
    cov2 <- EHF$cov2
    ## if (!is.null(cov1)) class(cov1) <- "matrix"
    ## if (!is.null(cov2)) class(cov2) <- "matrix"
    # }}}
    args <- list(object=object$crrFit,cov1=cov1,cov2=cov2,...)
    args <- args[!sapply(args,is.null)]
    ## warning("uuu")
    pred <- do.call(cmprsk::predict.crr,args)
    ## warning("bbb")
    out <- pred[,-1,drop=FALSE]
    if (!missing(times)){
        tind <- prodlim::sindex(jump.times=pred[,1],eval.times=times)
        out <- rbind(0,out)[tind+1,,drop=FALSE]
    }
    t(out)
}
