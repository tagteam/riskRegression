predict.FGR <- function(object,newdata,times,...){
  # {{{ check the data and the design
  if (missing(newdata)) stop("Argument 'newdata' is missing")
  formList <- readFormula(object$call$formula,specials=c("cov2"),specialArgumentNames=list("cov2"="tf"),unspecified="cov1")
  if (NROW(newdata) == 0) stop("No (non-missing) observations")
  # }}}
  # {{{ covariate design matrices
  cov1 <- modelMatrix(formula=formList$cov1$formula,
                   data=newdata,
                   intercept=NULL)
  cov2 <- modelMatrix(formula=formList$cov2$formula,
                      data=newdata,
                      intercept=NULL)
  if (!is.null(cov1))
    class(cov1) <- "matrix"
  if (!is.null(cov2))
    class(cov2) <- "matrix"
  # }}}
  args <- list(object=object$crrFit,cov1=cov1,cov2=cov2,...)
  args <- args[!sapply(args,is.null)]
  predcrr <- cmprsk:::predict.crr
  pred <- do.call("predcrr",args)
  out <- pred[,-1,drop=FALSE]
  if (!missing(times)){
    tind <- sindex(jump.times=pred[,1],eval.times=times)
    out <- rbind(0,out)[tind+1,,drop=FALSE]
  }
  t(out)
}
