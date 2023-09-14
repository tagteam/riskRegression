#' @title Fitting GLMnet for use with predictRisk
#'
#' @description Fit GLMnet models via a formula and a data set for use with \code{\link{predictRisk}}.
#' @name GLMnet
#'
#' @param formula A formula.
#' @param data The data on which to fit the model.
#' @param lambda The tuning parameters for GLMnet. If set to NULL, then it the parameters are chosen for you.
#' @param cv Whether to use cross-validation or not. Default is TRUE.
#' @param alpha The elasticnet mixing parameter. See the ?glmnet for more details.
#' @param nfolds Number of folds for cross-validation. Default is 10.
#' @param type.measure loss to use for cross-validation. Default is deviance.
#' @param family passed to \code{glmnet}. Defaults for binary outcome to \code{"binomial"} and for survival to \code{"cox"}.
#' @param \dots Additional arguments that are passed on to the glmnet.
#' @export
GLMnet <- function(formula,
                   data,
                   lambda = NULL,
                   cv = TRUE,
                   alpha = 1,
                   nfolds = 10,
                   type.measure = "deviance",
                   family,
                   ...) {
  requireNamespace(c("glmnet", "prodlim"))
  tt <- all.vars(update(formula, ".~1"))
  if (length(tt) != 1) {
    if (missing(family)) family <- "cox"
    strata.num <- start <- status <- NULL
    EHF <- prodlim::EventHistory.frame(formula, data, unspecialsDesign = TRUE, specials = NULL)
    stopifnot(attr(EHF$event.history, "model")[[1]] == "survival")
    # blank Cox object needed for predictions
    data <- data.frame(cbind(EHF$event.history, EHF$design))
    bl_cph <- coxph(Surv(time, status) ~ 1, data = data, x = 1, y = 1)
    bl_obj <- coxModelFrame(bl_cph)[]
    bl_obj[, strata.num := 0]
    bl_obj[, order := 1:.N]
    data.table::setorder(bl_obj, strata.num, stop, start, -status)
    sorted_x_train <- EHF$design[bl_obj[, order], ] ## save this for calculating Breslow estimator when doing predictions
    ## bl_obj$linear_predictor = c(exp(predict(fit,newx = EHF$design,type = "link", s=lambda)))
    design <- EHF$design[bl_obj$order, ]
    if (!cv) {
      fit <- glmnet::glmnet(
        x = EHF$design,
        y = EHF$event.history,
        lambda = lambda,
        alpha = alpha,
        family = family,
        ...
      )
    } else {
      fit <- glmnet::cv.glmnet(
        x = EHF$design,
        y = EHF$event.history,
        lambda = lambda,
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha,
        family = family,
        ...
      )
      lambda <- fit$lambda
    }
  } else {
    if (missing(family)) family <- "binomial"
    sorted_x_train <- bl_obj <- terms <- design <- NULL
    y <- data[[tt[1]]]
    x <- model.matrix(formula, data = data)
    if (!cv) {
      fit <- glmnet::glmnet(
        x = x,
        y = y,
        lambda = lambda,
        alpha = alpha,
        family = family,
        ...
      )
    } else {
      fit <- glmnet::cv.glmnet(
        x = x,
        y = y,
        lambda = lambda,
        nfolds = nfolds,
        type.measure = type.measure,
        alpha = alpha,
        family = family,
        ...
      )
      lambda <- fit$lambda
    }
  }
  out <- list(
    fit = fit,
    surv_info = bl_obj,
    design = design,
    call = match.call(),
    terms = terms(formula),
    sorted_x_train = sorted_x_train,
    cv = cv,
    lambda = lambda
  )
  class(out) <- "GLMnet"
  out
}
