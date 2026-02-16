#' @title Formula interface for glmnet 
#'
#' @description Fit glmnet models via a formula and a data set for use with \code{\link{predictRisk}}.
#' @name GLMnet
#'
#' @param formula Formula where the left hand side specifies either a single variable (continuous, binary or categorical),
#' or as a survival outcome (time, event), and the right hand side specifies the linear predictor.
#' Survival outcome can either be specified via \link[survival]{Surv} or via \link[prodlim]{Hist}.
#' Variables on the right hand side of the formula can be marked as \code{unpenalized}, see examples.
#' #' For survival outcome, a penalized Cox regression model is fitted. For this the formula may specify variables for which the baseline hazard function should be stratified, see examples. 
#' @param data The data used to fit the model.
#' @param lambda A hyperparameter passed to glmnet. If set to NULL, then lambda is chosen by cross-validation,
#' via the function \link[glmnet]{cv.glmnet}
#' @param alpha The elasticnet mixing parameter, with 0<=alpha<= 1. \code{alpha =1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#' @param cv Whether to use cross-validation or not. Default is TRUE.
#' @param nfolds Passed on to \link[glmnet]{cv.glmnet}. Number of folds for cross-validation. The default is 10.
#' @param type.measure Passed on to \link[glmnet]{cv.glmnet}. Loss to use for cross-validation. Default is deviance.
#' @param selector On of \code{'min'}, \code{'1se'}, \code{'undersmooth'} where the first two are described in the help page of \link[glmnet]{cv.glmnet}
#'                      and the latter is the smallest lambda value where the model could fit. Default is \code{'min'}.
#'        When \code{'undersmooth'} is specified no cross-validation is performed. 
#' @param family Passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}. For binary outcome the default is \code{"binomial"} and for survival \code{"cox"}.
#' @param \dots Additional arguments that are passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}.
#' @return A glmnet object enhanced with the call, the terms to create the design matrix, and in the survival case
#' with the Breslow estimate of the baseline hazard function. 
#' @seealso \link[glmnet]{glmnet}
#' @examples
#' library(survival)
#' set.seed(8)
#' d <- sampleData(87,outcome="survival")
#' test <- sampleData(5,outcome="survival")
#'
#' # penalized logistic regression
#' g <- GLMnet(X2~X1+X8,data=d,family="binomial")
#' predictRisk(g,newdata=test)
#' \dontrun{
#' g1 <- GLMnet(X2~X1+X8,data=d,lambda=0,gamma=0.5)
#' g2 <- GLMnet(X2~X1+X8,data=d,relax=TRUE)
#' }
#' # penalized Cox regression
#' 
#' f0 <- GLMnet(Surv(time,event)~X1+X2+X8+X9,data=d,lambda=0)
#' predictRisk(f0,newdata=test,times=3)
#' f <- GLMnet(Surv(time,event)~X1+X2+X8+X9,data=d)
#' f
#' predictCox(f,newdata=test,times=5,product.limit=TRUE)
#' predictRisk(f,newdata=test,times=1)
#' 
#' f1 <- GLMnet(Surv(time,event)~X1+X2+unpenalized(X8)+X9,data=d)
#' 
#' predictRisk(f1,newdata=test,times=1)
#' @export
GLMnet <- function(formula,
                   data,
                   lambda=NULL,
                   alpha = 1,
                   cv=TRUE,
                   nfolds = 10,
                   type.measure = "deviance",
                   selector = "min",
                   family,
                   ...){
    requireNamespace(c("glmnet","prodlim"))
    response_variable_names <- all.vars(update(formula,".~1"))
    selector <- match.arg(selector,c("undersmooth","min","1se"))
    # disable cross-validation when there is nothing to choose from or
    # when the choice does not require a criterion
    if (length(lambda) == 1 || selector == "undersmooth") cv <- FALSE
    if (length(lambda) == 1) selector <- "prespecified"
    if (length(response_variable_names) != 1){
        response <- prodlim::EventHistory.frame(formula = formula,
                                                data = data,
                                                unspecialsDesign = TRUE,
                                                specialsDesign = TRUE,
                                                stripSpecials = c("unpenalized","strata"),
                                                specials = c("strata","unpenalized"))
        Y <- response$event.history
        if (missing(family)) family <- "cox"
        if (attr(Y,"model")[[1]] != "survival"){
            stop("This function works only for survival models without competing risks.\nFor competing risks use riskRegression::CSC.")
        }
        if(length(response$strata)>0){
            # FIXME: can strata be based on multiple variables?
            attr(Y,"strata") <- response$strata
            class(Y) <- c("stratifySurv",class(Y))
        }
    }else{
        response <- Publish::specialFrame(formula = formula,
                                          data = data,
                                          specials.design = TRUE,
                                          unspecials.design = TRUE,
                                          strip.specials = "unpenalized",
                                          specials = "unpenalized")
        Y <- unlist(response$response)
        if (missing(family)){
            if (length(unique(Y)) == 2){
                family <- "binomial"
            } else{
                family <- "gaussian"
            }
        }
    }
    glmnet_args <- list(...)
    if ("penalty.factor" %in% names(glmnet_args)){
        penalty.factor <- glmnet_args$penalty
        X <- response$design
    }else{
        if (!is.null(response$unpenalized)){
            penalty.factor <- c(rep(1,NCOL(response$design)),rep(0,NCOL(response$unpenalized)))
            X <- cbind(response$design,response$unpenalized)
        }else{
            X <- response$design
            penalty.factor <- rep(1,NCOL(response$design))
        }
    }
    if (!cv){
        fit <- glmnet::glmnet(x=X,
                              y=Y,
                              lambda=lambda,
                              alpha= alpha,
                              penalty.factor = penalty.factor,
                              family=family,
                              ...)
        fit$call$alpha <- alpha
        fit$call$penalty.factor <- penalty.factor
        fit$call$family <- family
        fit$call$lambda <- lambda
        if (selector == "prespecified"){
            selected.lambda <- lambda
            selected.beta <- fit$beta
        }else{
            # undersmoothing
            selected.lambda <- fit$lambda[length(fit$lambda)]
            selected.beta <- fit$beta[,match(selected.lambda,fit$lambda,nomatch = NCOL(fit$beta)),drop = FALSE]
        }
    } else {
        # forcing cv
        fit <- glmnet::cv.glmnet(x=X,
                                 y=Y,
                                 lambda=lambda,
                                 nfolds=nfolds,
                                 type.measure = type.measure,
                                 alpha=alpha,
                                 penalty.factor = penalty.factor,
                                 family=family,
                                 ...)
        fit$call$alpha <- alpha
        fit$call$penalty.factor <- penalty.factor
        fit$call$type.measure <- type.measure
        fit$call$family <- family
        fit$call$nfolds <- nfolds
        fit$call$lambda <- lambda
        selected.lambda <- fit[[paste0("lambda.",selector)]]
        selected.beta <- fit$glmnet.fit$beta[,fit$index[,"Lambda"][[selector]],drop = FALSE]
    }
    colnames(selected.beta) <- switch(family,"binomial" = "log_odds",
                                      "cox" = "log_hazard_ratio",
                                      "beta")
    out <- list(fit = fit,
                y = Y,
                x = X,
                call = match.call(),
                terms = terms(formula,specials = c("strata","unpenalized")),
                cv=cv,
                selector = selector,
                penalty.factor = penalty.factor,
                selected.lambda = selected.lambda,
                selected.beta = selected.beta)
    class(out) <- "GLMnet"
    out
}
