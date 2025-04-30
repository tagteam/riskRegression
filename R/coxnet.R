### coxnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (08:29) 
## Version: 
## Last-Updated: Apr 29 2025 (14:32) 
##           By: Thomas Alexander Gerds
##     Update #: 7
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
#' history and the right hand side the linear predictor.  
##' @param data The data on which to fit the model. 
##' @param lambda The tuning parameters for GLMnet. If set to NULL, then it the parameters are chosen for you.
##' @param alpha
##' @param cv
##' @param strata
##' @param ... 
##' @return 
##' @seealso 
##' @examples 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
coxnet <- function(formula,
                   data,
                   lambda=0.1,
                   alpha=1,
                   cv=FALSE,
                   strata=NULL,
                   ...){
    # alpha = 1 equals lasso, alpha = 0 equals ridge, and alpha in (0,1) equals elastic net
    ehf <- prodlim::EventHistory.frame(formula, data, specials = NULL)
    if(length(strata)>0){
        ehf$event.history <- glmnet::stratifySurv(ehf$event.history, strata = data[[strata]])
    }
    if(cv){
        gfit <- glmnet::cv.glmnet(y = ehf$event.history, x = ehf$design, alpha = alpha, family = "cox",...)
    } else{
        gfit <- glmnet::glmnet(y = ehf$event.history, x = ehf$design, alpha = alpha, lambda = lambda, family = "cox",...)
    }
    gfit$design <- ehf$design
    gfit$y <- ehf$event.history
    gfit$formula <- formula
    gfit$strata <- strata
    gfit$cv <- cv
    class(gfit) <- c(class(gfit), "coxnet")
    gfit
}

######################################################################
### coxnet.R ends here
