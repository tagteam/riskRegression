#' @S3method predict FGR
predict.FGR <- function(object,newdata,times,...){
    # {{{ check the data and the design
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (NROW(newdata) == 0) stop("No (non-missing) observations")
    rhs <- as.formula(delete.response(object$terms))
    EHF <- prodlim::EventHistory.frame(rhs,
                                       newdata,
                                       check.formula=FALSE,
                                       specials=c("cov1","cov2"),
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
    pred <- do.call("predict.crr",args)
    out <- pred[,-1,drop=FALSE]
    if (!missing(times)){
        tind <- prodlim::sindex(jump.times=pred[,1],eval.times=times)
        out <- rbind(0,out)[tind+1,,drop=FALSE]
    }
    t(out)
}
