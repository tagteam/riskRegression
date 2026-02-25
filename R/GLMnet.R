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
#' For \code{rcs()}, GLMnet treats the term as a marker and will compute the spline basis using
#' \code{rms::rcspline.eval} with knots derived from the training data (or explicit knots supplied as the
#' second argument). To ensure GLMnet controls the knot choice, \code{rcs()} should refer to the marker
#' function \code{riskRegression::rcs} (which returns its input unchanged), not \code{rms::rcs}.
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
#' @param \dots Additional arguments passed on to \link[glmnet]{glmnet} and \link[glmnet]{cv.glmnet}.
#'   To avoid clashes, if \code{alpha} or \code{penalty.factor} are supplied in \code{...}, they are removed from \code{...}
#'   and handled explicitly. The final penalty vector used by GLMnet (after applying \code{unpenalized()},
#'   \code{pen()}, and any user-supplied \code{penalty.factor}) is stored as \code{penalty.factor.used}.
#'
#' @return A GLMnet object (class \code{"GLMnet"}) containing:
#' \itemize{
#'   \item \code{fit}: glmnet or cv.glmnet fit
#'   \item \code{x}: design matrix used for fitting
#'   \item \code{y}: response used for fitting
#'   \item \code{terms}: terms object with specials
#'   \item \code{penalty.factor.used}: named vector aligned with \code{colnames(x)}
#'   \item \code{rcs.info}: list of spline specs (knots + column mapping) for transforming new data consistently
#' }
#'
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
                   ...) {

  requireNamespace(c("glmnet", "prodlim"))

  .parse_pen_calls <- function(formula) {
    tt <- terms(formula)
    labs <- attr(tt, "term.labels")
    res <- list()
    if (length(labs) == 0) return(res)

    for (lab in labs) {
      expr <- try(str2lang(lab), silent = TRUE)
      if (inherits(expr, "try-error")) next
      if (!is.call(expr)) next
      if (!identical(expr[[1]], as.name("pen"))) next
      if (length(expr) < 3) stop("pen(x,w): please supply a penalty weight w.")
      vname <- as.character(expr[[2]])
      w <- as.numeric(eval(expr[[3]], parent.frame()))
      if (!is.finite(w) || length(w) != 1) stop("pen(x,w): w must be a single finite number.")
      res[[vname]] <- w
    }
    res
  }

  .parse_rcs_calls <- function(formula) {
    tt <- terms(formula)
    labs <- attr(tt, "term.labels")
    res <- list()
    if (length(labs) == 0) return(res)

    for (lab in labs) {
      expr <- try(str2lang(lab), silent = TRUE)
      if (inherits(expr, "try-error")) next
      if (!is.call(expr)) next
      if (!identical(expr[[1]], as.name("rcs"))) next
      vname <- as.character(expr[[2]])
      res[[vname]] <- expr
    }
    res
  }

  .rcs_knots_from_call <- function(call_obj, x_train) {
    requireNamespace("rms")
    if (length(call_obj) >= 3) {
      arg3 <- eval(call_obj[[3]], parent.frame())
      if (is.numeric(arg3) && length(arg3) == 1) {
        nk <- as.integer(arg3)
        if (!is.finite(nk) || nk < 3) stop("rcs(x,nk): nk must be >= 3.")
        knots <- Hmisc::rcspline.eval(x_train, nk = nk, knots.only = TRUE)
        return(list(knots = knots, nk = nk))
      } else if (is.numeric(arg3) && length(arg3) >= 3) {
        return(list(knots = as.numeric(arg3), nk = length(arg3)))
      } else {
        stop("rcs(x, ...): second argument must be nk (single number) or a numeric vector of knots.")
      }
    } else {
      nk <- 5L
      knots <- Hmisc::rcspline.eval(x_train, nk = nk, knots.only = TRUE)
      return(list(knots = knots, nk = nk))
    }
  }

  .apply_pen_weights_to_columns <- function(pf, cn, pen_map) {
    if (length(pen_map) == 0) return(pf)
    for (vname in names(pen_map)) {
      w <- pen_map[[vname]]
      hit <- (cn == vname) |
        startsWith(cn, paste0(vname)) |
        startsWith(cn, paste0(vname, ":")) |
        startsWith(cn, paste0("`", vname, "`"))
      if (any(hit)) pf[hit] <- pf[hit] * w
    }
    pf
  }

  response_variable_names <- all.vars(update(formula, ".~1"))
  selector <- match.arg(selector, c("undersmooth", "min", "1se"))

  if (length(lambda) == 1 || selector == "undersmooth") cv <- FALSE
  if (length(lambda) == 1) selector <- "prespecified"

  glmnet_args <- list(...)

  # (a) avoid alpha / penalty.factor clashes: pull them out of ...
  if ("alpha" %in% names(glmnet_args)) {
    alpha <- glmnet_args$alpha
    glmnet_args$alpha <- NULL
  }
  user_penalty_factor <- NULL
  if ("penalty.factor" %in% names(glmnet_args)) {
    user_penalty_factor <- glmnet_args$penalty.factor
    glmnet_args$penalty.factor <- NULL
  }

  pen_map <- .parse_pen_calls(formula)
  rcs_map <- .parse_rcs_calls(formula)

  # build response + design parts; include rcs as special so we can grab raw x from response$rcs
  if (length(response_variable_names) != 1) {

    response <- prodlim::EventHistory.frame(
      formula = formula,
      data = data,
      unspecialsDesign = TRUE,
      specialsDesign = TRUE,
      stripSpecials = c("unpenalized", "strata", "pen", "rcs"),
      specials = c("strata", "unpenalized", "pen", "rcs")
    )
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

    response <- Publish::specialFrame(
      formula = formula,
      data = data,
      specials.design = TRUE,
      unspecials.design = TRUE,
      strip.specials = c("unpenalized", "pen", "rcs"),
      specials = c("unpenalized", "pen", "rcs")
    )
    Y <- unlist(response$response)

    if (missing(family)) {
      if (length(unique(Y)) == 2) family <- "binomial" else family <- "gaussian"
    }
  }

  # ---- build X and base penalty.factor from pieces ---------------------------

  # rcs raw columns (as specials.design) should be penalized by default
  rcs_raw <- NULL
  if (!is.null(response$rcs)) {
    rcs_raw <- response$rcs
    # make sure rcs specials are single-column per variable (marker rcs should do that)
    # and name columns by the underlying variable names if possible
    if (length(rcs_map) > 0 && NCOL(rcs_raw) == length(rcs_map)) {
      colnames(rcs_raw) <- names(rcs_map)
    }
  }

  if (!is.null(response$unpenalized)) {
    X <- cbind(response$design, rcs_raw, response$unpenalized)
    penalty.factor <- c(rep(1, NCOL(response$design)),
                        if (is.null(rcs_raw)) numeric(0) else rep(1, NCOL(rcs_raw)),
                        rep(0, NCOL(response$unpenalized)))
  } else {
    X <- cbind(response$design, rcs_raw)
    penalty.factor <- c(rep(1, NCOL(response$design)),
                        if (is.null(rcs_raw)) numeric(0) else rep(1, NCOL(rcs_raw)))
  }

  X <- as.matrix(X)
  colnames_X <- colnames(X)
  if (!is.null(colnames_X)) names(penalty.factor) <- colnames_X

  # incorporate user penalty.factor (from ...) safely (bug fixed: use penalty.factor, not $penalty)
  if (!is.null(user_penalty_factor)) {
    if (!is.null(names(user_penalty_factor))) {
      upf <- rep(1, NCOL(X))
      names(upf) <- colnames(X)
      upf[names(user_penalty_factor)] <- as.numeric(user_penalty_factor)
      user_penalty_factor <- upf
    }
    if (length(user_penalty_factor) != NCOL(X)) {
      stop("penalty.factor (from ...) must have length ncol(design matrix) or be a named vector matching colnames(x).")
    }
    penalty.factor <- penalty.factor * as.numeric(user_penalty_factor)
    names(penalty.factor) <- colnames(X)
  }

  # ---- (d) expand rcs() using training knots and store rcs.info --------------

  rcs_info <- list()
  if (length(rcs_map) > 0) {
    requireNamespace("rms")
    for (vname in names(rcs_map)) {
      if (!vname %in% colnames(X)) next  # raw column must be present (from response$rcs)
      x_train <- data[[vname]]
      if (!is.numeric(x_train)) stop(sprintf("rcs(%s,...): only numeric variables are supported.", vname))

      hit <- which(colnames(X) == vname)
      if (length(hit) != 1) stop(sprintf("rcs(%s,...): expected exactly one raw column in X.", vname))

      spec <- .rcs_knots_from_call(rcs_map[[vname]], x_train)
      knots <- spec$knots

      basis <- Hmisc::rcspline.eval(x_train, knots = knots, inclx = TRUE)
      bcn <- paste0(vname, "_rcs", seq_len(NCOL(basis)))
      colnames(basis) <- bcn

      # replace raw column by basis
      X <- cbind(X[, -hit, drop = FALSE], basis)
      raw_pf <- penalty.factor[hit]
      penalty.factor <- c(penalty.factor[-hit], rep(raw_pf, NCOL(basis)))
      colnames(X) <- c(colnames_X[-hit], bcn)
      names(penalty.factor) <- colnames(X)
      colnames_X <- colnames(X)

      rcs_info[[vname]] <- list(variable = vname, knots = knots, inclx = TRUE, columns = bcn)
    }
  }

  # (b) apply pen() weights (after rcs expansion so weights apply to basis columns too)
  penalty.factor <- .apply_pen_weights_to_columns(penalty.factor, colnames(X), pen_map)
  names(penalty.factor) <- colnames(X)

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
    penalty.factor.used = penalty.factor,
    rcs.info = rcs_info,
    selected.lambda = selected.lambda,
    selected.beta = selected.beta
  )

  class(out) <- "GLMnet"
  out
}
