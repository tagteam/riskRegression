#' @title Fitting HAL for use with predictRisk
#'
#' @description Fit HAL models via a formula and a data set for use with \code{\link{predictRisk}}.
#' @name Hal9001
#'
#' @param formula A formula.
#' @param data The data on which to fit the model. 
#' @param lambda The tuning parameters for HAL. If set to NULL, then it the parameters are chosen for you.
#' @param \dots Additional arguments that are passed on to the hal_fit.
#' @export
Hal9001 <- function(formula,data,lambda=NULL,...){
  requireNamespace("hal9001")
  strata.num = start = status = NULL
  EHF = prodlim::EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
  stopifnot(attr(EHF$event.history,"model")[[1]] == "survival")
  # blank Cox object needed for predictions
  data = data.frame(cbind(EHF$event.history,EHF$design))
  bl_cph <- coxph(Surv(time,status)~1,data=data,x=1,y=1)
  bl_obj <- coxModelFrame(bl_cph)[]
  bl_obj[,strata.num:=0]
  data.table::setorder(bl_obj, strata.num,stop,start,-status)
  hal_fit <- do.call(hal9001::fit_hal,list(X = EHF$design,
                                           Y = EHF$event.history,
                                           family = "cox",
                                           return_lasso = TRUE,
                                           yolo = FALSE,
                                           lambda = lambda,
                                           ...))
  out = list(fit = hal_fit,
             surv_info = bl_obj,
             call = match.call(),
             terms = terms(formula))
  class(out) = "Hal9001"
  out
}

#' @title Fitting GLMnet for use with predictRisk
#'
#' @description Fit GLMnet models via a formula and a data set for use with \code{\link{predictRisk}}.
#' @name GLMnet
#'
#' @param formula A formula.
#' @param data The data on which to fit the model. 
#' @param lambda The tuning parameters for GLMnet. If set to NULL, then it the parameters are chosen for you.
#' @param cv Whether to use cross-validation or not. Default is TRUE.
#' @param alpha The elasticnet mixing parameter - is not used by cross-validation. See the ?glmnet for more details.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param type.measure loss to use for cross-validation. Default is deviance.
#' @param \dots Additional arguments that are passed on to the glmnet.
#' @export
GLMnet <- function(formula,data,lambda=NULL, cv=TRUE,alpha = 1,nfolds = 10, type.measure = "deviance",...){
  requireNamespace(c("glmnet","prodlim"))
  tt <- all.vars(update(formula,".~1"))
  if (length(tt) != 1){
    strata.num = start = status = NULL
    EHF = prodlim::EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
    stopifnot(attr(EHF$event.history,"model")[[1]] == "survival")
    # blank Cox object needed for predictions
    data = data.frame(cbind(EHF$event.history,EHF$design))
    bl_cph <- coxph(Surv(time,status)~1,data=data,x=1,y=1)
    bl_obj <- coxModelFrame(bl_cph)[]
    bl_obj[,strata.num:=0]
    data.table::setorder(bl_obj, strata.num,stop,start,-status)
    if (!cv){
      fit <- glmnet::glmnet(x=EHF$design,y=EHF$event.history, lambda=lambda, alpha= alpha,family="cox")
    }
    else {
      fit <- glmnet::cv.glmnet(x=EHF$design,y=EHF$event.history, lambda=lambda,nfolds=nfolds,type.measure = "deviance", family="cox")
      lambda <- fit$lambda
    }
  }
  else {
    bl_obj=terms=NULL
    y  <- data[[tt[1]]]
    x <- model.matrix(formula, data=data)
    if (!cv){
      fit <- glmnet::glmnet(x=x,y=y, lambda=lambda,alpha=alpha,family="binomial")
    }
    else {
      fit <- glmnet::cv.glmnet(x=x,y=y, lambda=lambda,nfolds=nfolds,type.measure =type.measure, family="binomial")
      lambda <- fit$lambda
    }
  }
  out = list(fit = fit,
             surv_info = bl_obj,
             call = match.call(),
             terms = terms(formula),
             cv=cv,
             lambda = lambda)
  class(out) = "GLMnet"
  out
}
