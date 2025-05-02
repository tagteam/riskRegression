### coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (08:29) 
## Version: 
## Last-Updated: May  2 2025 (12:32) 
##           By: Thomas Alexander Gerds
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Formula interface for penalized Cox regression models via glmnet
##'
##' This is a wrapper function which fits a Cox regression model via glmnet
##' and enhances the resulting object 
##' @title
##' @param formula Formula where the left hand side specifies the event
##' history and the right hand side the linear predictor.  
##' @param data The data on which to fit the model. 
##' @param strata
##' @param lambda The tuning parameters for GLMnet. If set to NULL, then it the parameters are chosen for you.
#' @param alpha The elasticnet mixing parameter. 
##' @param nfolds Number of folds for cross-validation. Default is 10.
##' @param ...  Additional arguments that are passed on to glmnet.
##' @return cox net object
##' @seealso \link[glmnet]{glmnet}
##' @examples
##' d <- sampleData(73,outcome="survival")
##' f <- coxnet(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=d)
##' f1 <- coxnet(Surv(time,event)~strata(X1)+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=d)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
coxnet <- function(formula,
                   data,
                   lambda=NULL,
                   alpha=1,
                   cv=FALSE,
                   nfolds = 10,
                   ...){
    requireNamespace("glmnet")
    unpenalized <- function(x)x
    # alpha = 1 equals lasso, alpha = 0 equals ridge, and alpha in (0,1) equals elastic net
    EHF <- prodlim::EventHistory.frame(formula, data, specials = c("strata","unpenalized"))
    if(length(EHF$strata)>0){
        attr(EHF$event.history,"strata") <- EHF$strata
        class(EHF$event.history) <- c("stratifySurv",class(EHF$event.history))
    }
    if (length(EHF$unpenalized)>0){
        penalty.factor <- c(rep(1,NCOL(EHF$design)),rep(0,NCOL(EHF$unpenalized)))
        X <- cbind(EHF$design,EHF$unpenalized)
    }else{
        X <- EHF$design
        penalty.factor <- rep(1,NCOL(EHF$design))
    }
    args <- list(y = EHF$event.history,
                 x = X,
                 alpha = alpha,
                 family = "cox",
                 nfolds=nfolds,
                 penalty.factor = penalty.factor,
                 ...)
    if (cv){
        cv.glmnet(y = EHF$event.history,
                  x = X,
                  alpha = alpha,
                  family = "cox",
                  nfolds=nfolds,
                  penalty.factor = penalty.factor,
                  ...)
    }else{
        args$lambda <- lambda
        fun <- "glmnet"
    }
    fit <- do.call(what = fun ,args = args)
    fit$design <- EHF$design
    fit$strata <- strata
    fit$terms = terms(formula)
    fit
}

######################################################################
### coxnet.R ends here
