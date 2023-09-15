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
Hal9001 <- function(formula, data, lambda = NULL, ...) {
  requireNamespace("hal9001")
  strata.num <- start <- status <- NULL
  EHF <- prodlim::EventHistory.frame(formula, data, unspecialsDesign = TRUE, specials = NULL)
  stopifnot(attr(EHF$event.history, "model")[[1]] == "survival")
  # blank Cox object needed for predictions
  data <- data.frame(cbind(EHF$event.history, EHF$design))
  bl_cph <- coxph(Surv(time, status) ~ 1, data = data, x = 1, y = 1)
  bl_obj <- coxModelFrame(bl_cph)[]
  bl_obj[, strata.num := 0]
  data.table::setorder(bl_obj, strata.num, stop, start, -status)
  hal_fit <- do.call(hal9001::fit_hal, list(
    X = EHF$design,
    Y = EHF$event.history,
    family = "cox",
    return_lasso = TRUE,
    yolo = FALSE,
    lambda = lambda,
    ...
  ))
  out <- list(
    fit = hal_fit,
    surv_info = bl_obj,
    call = match.call(),
    terms = terms(formula)
  )
  class(out) <- "Hal9001"
  out
}
