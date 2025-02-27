#' Predict subject specific risks (cumulative incidence) based on Fine-Gray regression model
#'
#' Predict subject specific risks (cumulative incidence) based on Fine-Gray regression model
#' @param object Result of call to \code{FGR}
#' @param newdata Predictor values of subjects for who to predict risks
#' @param times Time points at which to evaluate the risks 
#' @param ... passed to predict.crr
#'
#' @examples
##' library(prodlim)
#' library(survival)
#' set.seed(10)
#' d <- sampleData(101, outcome = "competing.risk")
#' tFun<-function(t) {t}
#' fgr<-FGR(Hist(time, event)~X1+strata(X2)+X6+cov2(X7, tf=tFun),
#'          data=d, cause=1)
#' predictRisk(fgr,times=5,newdata=d[1:10])
#' @method predict FGR
#' @export
predict.FGR <- function(object,newdata,times,...){
    # {{{ check the data and the design
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (NROW(newdata) == 0) stop("No (non-missing) observations")
    rhs <- as.formula(delete.response(object$terms))
    ## EHF <- prodlim::model.design(rhs,
                                 ## newdata,
                                 ## dropIntercept=TRUE,
    ## specialsDesign=TRUE)
    ## need "response" to make EventHistory.frame work
    ## but Score removes response for prediction
    newdata$dummy.time=1
    newdata$dummy.event=1
    dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    EHF <- prodlim::EventHistory.frame(dummy.formula,
                                       data=newdata,
                                       specials=c("cov1","cov2"),
                                       stripSpecials=c("cov1","cov2"),
                                       stripArguments=list("cov1"=NULL,"cov2"=list("tf"=NULL)),
                                       stripUnspecials="cov1",
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
    ## print(colnames(args$cov1))
    pred <- do.call(cmprsk::predict.crr,args)
    ## warning("bbb")
    out <- pred[,-1,drop=FALSE]
    if (!missing(times)){
        tind <- prodlim::sindex(jump.times=pred[,1],eval.times=times)
        out <- rbind(0,out)[tind+1,,drop=FALSE]
    }
    t(out)
}
