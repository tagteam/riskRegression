#' @title Fitting HAL for use with predictRisk
#'
#' @description Fit HAL models via a formula and a data set for use with \code{\link{predictRisk}}.
#' @name Hal9001
#' @param formula A formula.
#' @param data The data on which to fit the model. 
#' @param lambda The tuning parameters for HAL. If set to NULL, then it the parameters are chosen for you.
#' @param family Can be one of \code{c("binomial","normal")}. If any other value is translated into event history analysis.
#' @param \dots Additional arguments that are passed on to the hal_fit.
#' @export
Hal9001 <- function(formula,data,lambda=NULL,family,...){
    requireNamespace("hal9001")
    strata.num = start = status = NULL
    if (family %in% c("binomial","normal")){
        response <- Publish::specialFrame(formula = formula,
                                          data = data,
                                          specials.design = TRUE,
                                          unspecials.design = TRUE,
                                          strip.specials = "unpenalized",
                                          specials = "unpenalized")
        Y <- unlist(response$response)
        X <- cbind(response$design,response$unpenalized)
        hal_fit <- do.call(hal9001::fit_hal,list(X = X,
                                                 Y = Y,
                                                 family = family,
                                                 return_lasso = TRUE,
                                                 yolo = FALSE,
                                                 lambda = lambda,
                                                 ...))
        out = list(fit = hal_fit,
                   call = match.call(),
                   terms = terms(formula))
    }else{ 
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
        # to be able to calculate the baseline hazard function we need
        # to evaluate the linear predictor with the training covariates
        data.table::set(bl_obj,j = "eXb",value= as.numeric(predict(hal_fit,new_data=EHF$design)))
        out = list(fit = hal_fit,
                   surv_info = bl_obj,
                   call = match.call(),
                   terms = terms(formula))
    }
    class(out) = "Hal9001"
    out
}

