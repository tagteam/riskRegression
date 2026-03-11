#' @title Formula interface for glmnet
#'
#' @description
#' Fit glmnet models via a formula and a data set for use with \code{\link{predictRisk}}.
#'
#' Supports special terms on the right-hand side:
#' \itemize{
#'   \item \code{unpenalized(x)}: exclude variable \code{x} from penalization.
#'   \item \code{pen(x, w)}: apply a relative penalty weight \code{w} to variable \code{x}
#'         (e.g., \code{pen(x,2)} more penalty; \code{pen(z,0)} unpenalized).
#'   \item \code{rcs(x, ...)}: restricted cubic splines for numeric \code{x}. This is handled as a
#'         *special marker* (see Details). Knots are estimated from the training data (unless specified)
#'         and stored in the fitted object.
#'   \item \code{strata(x)} (survival outcome only): stratified baseline hazards.
#' }
#'
#' @details
#' For \code{rcs()}, GLMnet will compute the spline basis using
#' \link[rms]{rcs} with knots derived from the data. In this application, rcs has
#' two arguments: the number of knots \code{nknots} (defaults to 3) and the penalty factor \code{pf} (defaults to 2).
#'
#' @param formula Formula where the left hand side specifies either a single variable (continuous, binary or categorical),
#'   or a survival outcome (time, event), and the right hand side specifies the linear predictor.
#'   Survival outcome can either be specified via \link[survival]{Surv} or via \link[prodlim]{Hist}.
#' @param data The data used to fit the model.
#' @param lambda A hyperparameter passed to glmnet. If set to NULL, then lambda is chosen by cross-validation,
#'   via the function \link[glmnet]{cv.glmnet}.
#' @param alpha The elasticnet mixing parameter, with 0<=alpha<= 1. \code{alpha =1} is the lasso penalty, and \code{alpha=0} the ridge penalty.
#'   If \code{alpha} is also provided via \code{...}, the value in \code{...} takes precedence.
#' @param cv Whether to use cross-validation or not. Default is TRUE.
#' @param nfolds Passed on to \link[glmnet]{cv.glmnet}. Number of folds for cross-validation. The default is 10.
#' @param type.measure Passed on to \link[glmnet]{cv.glmnet}. Loss to use for cross-validation. Default is deviance.
#' @param selector One of \code{'min'}, \code{'1se'}, \code{'undersmooth'} where the first two are described in the help page of \link[glmnet]{cv.glmnet}
#'   and the latter is the smallest lambda value where the model could fit. Default is \code{'min'}.
#'   When \code{'undersmooth'} is specified, no cross-validation is performed.
#' @param family Passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}. For binary outcome the default is \code{"binomial"} and for survival \code{"cox"}.
#' @param penalty.factor The penalty factor which is passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}. It must be a named
#' vector where the names are the column names of the design matrix. To see the column names of the design matrix perform a dryrun without setting
#' argument \code{penalty.factor} and check element \code{fit$penalty.factor}.
#' @param \dots Additional arguments passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}.
#'
#' @return A GLMnet object (class \code{"GLMnet"}) containing:
#' \itemize{
#'   \item \code{fit}: glmnet or cv.glmnet fit
#'   \item \code{x}: design matrix used for fitting
#'   \item \code{y}: response used for fitting
#'   \item \code{terms}: terms object with specials
#'   \item \code{penalty.factor}: named vector aligned with \code{colnames(x)}
#'   \item \code{rcs.parameters}: named vector with parameters for the rcs terms
#' }
#' @examples
#' library(glmnet)
#' library(riskRegression)
#' library(survival)
#' data(Melanoma)
#' fit <- GLMnet(Surv(time,status==1)~pen(invasion,4)+epicel+ulcer+logthick+pen(sex,0)+age,
#'              data=Melanoma)
#' @seealso \link[glmnet]{glmnet}
#' @export
GLMnet <- function(formula,
                   data,
                   lambda = NULL,
                   alpha = 1,
                   cv = TRUE,
                   nfolds = 10,
                   type.measure = "deviance",
                   selector = "min",
                   family,
                   penalty.factor,
                   ...) {
    requireNamespace(c("glmnet", "prodlim"))
    response_variable_names <- all.vars(update(formula, ".~1"))
    selector <- match.arg(selector, c("undersmooth", "min", "1se"))
    # 
    if (length(lambda) == 1 || selector == "undersmooth") cv <- FALSE
    if (length(lambda) == 1) selector <- "prespecified"
    # 
    glmnet_args <- list(...)
    # 
    # alpha / penalty.factor clashes: pull them out of ...
    if ("alpha" %in% names(glmnet_args)) {
        alpha <- glmnet_args$alpha
        glmnet_args$alpha <- NULL
    }
    if ("penalty.factor" %in% names(glmnet_args)) {
        glmnet_args$penalty.factor <- NULL
    }
    # build response + design parts
    response <- prodlim::EventHistory.frame(
                             formula = formula,
                             # do not check for survival when Y is a single variable 
                             check.formula = length(response_variable_names) > 1,
                             data = data,
                             unspecialsDesign = TRUE,
                             specialsDesign = TRUE,
                             stripSpecials = c("unpenalized", "strata", "pen", "rcs"),
                             stripArguments=list("pen"=list("pf" = 1),"unpenalized"=NULL,"rcs"=list("nknots" = 3,"pen" = 2),"strata" = NULL),
                             specials = c("strata", "unpenalized", "pen", "rcs"))
    if (length(response_variable_names) > 1) {
        Y <- response$event.history
        if (missing(family)) family <- "cox"
        if (attr(Y, "model")[[1]] != "survival") {
            stop("This function works only for survival models without competing risks.\nFor competing risks use riskRegression::CSC.")
        }
        if (length(response$strata) > 0) {
            attr(Y, "strata") <- response$strata
            class(Y) <- c("stratifySurv", class(Y))
        }
    } else {
        Y <- unlist(response$response)
        if (missing(family)) {
            if (length(unique(Y)) == 2) family <- "binomial" else family <- "gaussian"
        }
    }
    #
    # Design matrix
    #
    # build design matrix X and penalty.factor
    X <- cbind(response$design, response$pen, response$unpenalized)
    colnames_X <- colnames(X)
    internal_penalty.factor <- sapply(colnames_X,function(x){1})
    rcs.parameters <- NULL
    # add restricted cubic spline terms
    if (!is.null(response$rcs)) {
        rcs_knots <- attr(response$rcs,"arguments.terms")$nknots
        rcs_pen <- attr(response$rcs,"arguments.terms")$pen
        for (v in names(rcs_knots)){
            v_X <- do.call(rms::rcs,list(x = data[[v]],parms = as.numeric(rcs_knots[[v]])))
            v_parms <- list(attr(v_X,"parms"))
            names(v_parms) <- v
            rcs.parameters <- c(rcs.parameters,v_parms)
            colnames(v_X) <- paste0("rcs_",v,"_term",1:NCOL(v_X))
            v_penalty.factor <- sapply(colnames(v_X),function(u){rcs_pen[[v]]})
            internal_penalty.factor <- c(internal_penalty.factor,
                                         v_penalty.factor)
            X <- cbind(X,v_X)
        }
    }
    if (!is.null(response$unpenalized)) {
        unpenalized_terms <- unlist(attr(response$unpenalized,"matrix.terms"))
        internal_penalty.factor[intersect(unpenalized_terms,names(internal_penalty.factor))] <- 0
    }
    if (!is.null(response$pen)) {
        # pen uses argument pf which defaults to 1 as specified above
        formula_penalty_terms <- sapply(attr(response$pen,"arguments.terms")$pf,function(x){as.numeric(x)})
        for (this_term in intersect(names(internal_penalty.factor),names(formula_penalty_terms))){
            if (is.na(formula_penalty_terms[[this_term]]))stop(paste0("Penalty factor of term ",this_term," cannot be coerced to numeric"))
            internal_penalty.factor[this_term] <- formula_penalty_terms[[this_term]]
        }
    }
    if (missing(penalty.factor)){
        penalty.factor <- internal_penalty.factor
    }else{
        ## check user supplied penalty factor
        if (!(all(names(penalty.factor)%in% colnames_X))){
            stop(paste0("riskRegression::GLMnet: Argument penalty.factor must be a named vector which matches the column names of the design matrix:\n",
                        paste0(colnames_X,collapse = ", ")))
        }
        internal_penalty.factor[intersect(names(internal_penalty.factor),names(penalty.factor))] <- penalty.factor[intersect(names(internal_penalty.factor),names(penalty.factor))]
        penalty.factor <- internal_penalty.factor
    }
    # ---- fit glmnet / cv.glmnet without argument clashes -----------------------
    if (!cv) {
        fit <- do.call(glmnet::glmnet, c(list(
                                           x = X,
                                           y = Y,
                                           lambda = lambda,
                                           alpha = alpha,
                                           penalty.factor = penalty.factor,
                                           family = family
                                       ), glmnet_args))

        if (selector == "prespecified") {
            selected.lambda <- lambda
            selected.beta <- fit$beta
        } else {
            selected.lambda <- fit$lambda[length(fit$lambda)]
            selected.beta <- fit$beta[, match(selected.lambda, fit$lambda, nomatch = NCOL(fit$beta)), drop = FALSE]
        }

    } else {
        fit <- do.call(glmnet::cv.glmnet, c(list(
                                              x = X,
                                              y = Y,
                                              lambda = lambda,
                                              nfolds = nfolds,
                                              type.measure = type.measure,
                                              alpha = alpha,
                                              penalty.factor = penalty.factor,
                                              family = family
                                          ), glmnet_args))

    selected.lambda <- fit[[paste0("lambda.", selector)]]
    selected.beta <- fit$glmnet.fit$beta[, fit$index[, "Lambda"][[selector]], drop = FALSE]
  }
    colnames(selected.beta) <- switch(family,
                                      "binomial" = "log_odds",
                                      "cox" = "log_hazard_ratio",
                                      "beta")
    out <- list(
        fit = fit,
        y = Y,
        x = X,
        call = match.call(),
        terms = terms(formula, specials = c("strata", "unpenalized", "pen", "rcs")),
        cv = cv,
        selector = selector,
        penalty.factor = penalty.factor,
        rcs.parameters = rcs.parameters,
        selected.lambda = selected.lambda,
        selected.beta = selected.beta
    )
  class(out) <- "GLMnet"
  out
}
