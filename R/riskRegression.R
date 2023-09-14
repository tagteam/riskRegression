#' Risk Regression
#' Fits a regression model for the risk of an event -- allowing for competing
#' risks.
#'
#' This is a wrapper for the function \code{comp.risk} from the timereg package.
#' The main difference is one marks variables in the formula that should have a
#' time-dependent effect whereas in \code{comp.risk} one marks variables that
#' should have a time constant (proportional) effect.
#' @aliases riskRegression ARR LRR
#' @param formula Formula where the left hand side specifies the event
#' history event.history and the right hand side the linear predictor.  See
#' examples.
#' @param data The data for fitting the model in which includes all
#' the variables included in formula.
#' @param times Vector of times. For each time point in \code{times}
#' estimate the baseline risk and the timevarying coefficients.
#' @param link \code{"relative"} for the absolute risk regression
#' model.  \code{"logistic"} for the logistic risk regression model.
#' \code{"prop"} for the Fine-Gray regression model.
#' @param cause The cause of interest.
#' @param conf.int If \code{TRUE} return the iid decomposition, that
#' can be used to construct confidence bands for predictions.
#' @param cens.model Specified the model for the (conditional)
#' censoring distribution used for deriving weights (IFPW). Defaults
#' to "KM" (the Kaplan-Meier method ignoring covariates) alternatively
#' it may be "Cox" (Cox regression).
#' @param cens.formula Right hand side of the formula used for fitting
#' the censoring model.  If not specified the right hand side of
#' \code{formula} is used.
#' @param max.iter Maximal number of iterations.
#' @param conservative If \code{TRUE} use variance formula that ignores the contribution
#' by the estimate of the inverse of the probability of censoring weights
#' @param ... Further arguments passed to \code{comp.risk}
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}, Thomas H. Scheike \email{ts@@biostat.ku.dk}
#' @references
#' Thomas A Gerds, Thomas H Scheike, and Per K Andersen. Absolute risk
#' regression for competing risks: interpretation, link functions, and
#' prediction. Statistics in medicine, 31(29):3921--3930, 2012.
#'
#' Scheike, Zhang and Gerds (2008), Predicting cumulative incidence probability
#' by direct binomial regression, Biometrika, 95, 205-220.
#'
#' Scheike and Zhang (2007), Flexible competing risks regression modelling and
#' goodness of fit, LIDA, 14, 464-483.
#'
#' Martinussen and Scheike (2006), Dynamic regression models for survival data,
#' Springer.
##' @examples
##'
##' library(prodlim)
##' data(Melanoma,package="riskRegression")
##' ## tumor thickness on the log-scale
##' Melanoma$logthick <- log(Melanoma$thick)
##'
##' # Single binary factor
##'
##' ## absolute risk regression
##' library(survival)
##' library(prodlim)
##' fit.arr <- ARR(Hist(time,status)~sex,data=Melanoma,cause=1)
##' print(fit.arr)
##' # show predicted cumulative incidences
##' plot(fit.arr,col=3:4,newdata=data.frame(sex=c("Female","Male")))
##'
##' ## compare with non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit.aj <- prodlim(Hist(time,status)~sex,data=Melanoma)
##' plot(fit.aj,conf.int=FALSE)
##' plot(fit.arr,add=TRUE,col=3:4,newdata=data.frame(sex=c("Female","Male")))
##'
##' ## with time-dependent effect
##' fit.tarr <- ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1)
##' plot(fit.tarr,newdata=data.frame(sex=c("Female","Male")))
##'
##' ## logistic risk regression
##' fit.lrr <- LRR(Hist(time,status)~sex,data=Melanoma,cause=1)
##' summary(fit.lrr)
##'
##'
##' # Single continuous factor
##'
##' ## tumor thickness on the log-scale
##' Melanoma$logthick <- log(Melanoma$thick)
##'
##' ## absolute risk regression
##' fit2.arr <- ARR(Hist(time,status)~logthick,data=Melanoma,cause=1)
##' print(fit2.arr)
##' # show predicted cumulative incidences
##' plot(fit2.arr,col=1:5,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##'
##' ## comparison with nearest neighbor non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit2.aj <- prodlim(Hist(time,status)~logthick,data=Melanoma)
##' plot(fit2.aj,conf.int=FALSE,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##' plot(fit2.arr,add=TRUE,col=1:5,lty=3,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##'
##' ## logistic risk regression
##' fit2.lrr <- LRR(Hist(time,status)~logthick,data=Melanoma,cause=1)
##' summary(fit2.lrr)
##'
##' ## change model for censoring weights
##' library(rms)
##' fit2a.lrr <- LRR(Hist(time,status)~logthick,
##'                  data=Melanoma,
##'                  cause=1,
##'                  cens.model="cox",
##'                  cens.formula=~sex+epicel+ulcer+age+logthick)
##' summary(fit2a.lrr)
##'
##' ##  compare prediction performance
##' Score(list(ARR=fit2.arr,AJ=fit2.aj,LRR=fit2.lrr),formula=Hist(time,status)~1,data=Melanoma)
##'
##'
##' # multiple regression
##' library(riskRegression)
##' library(prodlim)
##' # absolute risk model
##' multi.arr <- ARR(Hist(time,status)~logthick+sex+age+ulcer,data=Melanoma,cause=1)
##'
##' # stratified model allowing different baseline risk for the two gender
##' multi.arr <- ARR(Hist(time,status)~thick+strata(sex)+age+ulcer,data=Melanoma,cause=1)
##'
##' # stratify by a continuous variable: strata(age)
##' multi.arr <- ARR(Hist(time,status)~tp(thick,power=0)+strata(age)+sex+ulcer,
##'                  data=Melanoma,
##'                  cause=1)
##'
##' fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
##' summary(fit.arr2a)
##' fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
##' summary(fit.arr2b)
##'
##' ## logistic risk model
##' fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
##' summary(fit.lrr)
##'
##'
##'
##'
##'
##' ## nearest neighbor non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit.aj <- prodlim(Hist(time,status)~thick,data=Melanoma)
##' plot(fit.aj,conf.int=FALSE)
##'
##' # prediction performance
##' x <- Score(list(fit.arr2a,fit.arr2b,fit.lrr),
##'              data=Melanoma,
##'              formula=Hist(time,status)~1,
##'              cause=1,
##'              split.method="none")
##'
##'
#' @keywords survival
#' @export
riskRegression <- function(formula,
                           data,
                           times,
                           link = "relative",
                           cause,
                           conf.int = TRUE,
                           cens.model,
                           cens.formula,
                           max.iter = 50,
                           conservative = TRUE,
                           ...) {
  # {{{ preliminaries
# weighted <- 0
# detail <- 0
  stopifnot(is.numeric(max.iter) && max.iter[[1]] > 0 && (round(max.iter[[1]]) == max.iter[[1]]))
  # trans=1 P_1=1-exp(- ( x' b(b)+ z' gam t^time.pow) ),
  # trans=2 P_1=1-exp(-exp(x a(t)+ z` b )
  # trans=not done P_1=1-exp(-x a(t) exp(z` b )) is not good numerically
  # trans=3 P_1=exp(-exp(x a(t)+ z` b )
  ##   trans <- switch(link,"log"=
  trans <- switch(link,
    "additive" = "additive", #
    "prop" = "prop", # Proportional hazards (Cox, FG)
    "logistic" = "logistic2", # Logistic absolute risks
    "relative" = "rcif2"
  ) # Relative absolute risks
  # }}}
  # {{{ read the data and the design
  call <- match.call()
  EHF <- prodlim::EventHistory.frame(formula,
    data,
    specials = c("timevar", "strata", "prop", "const", "tp"),
    stripSpecials = c("timevar", "prop"),
    stripArguments = list("prop" = list("power" = 0), "timevar" = list("test" = 0)),
    stripAlias = list("timevar" = c("strata"), "prop" = c("tp", "const")),
    stripUnspecials = "prop",
    specialsDesign = TRUE,
    dropIntercept = TRUE
  )
  ## Terms <- attr(EHF,"Terms")
  Terms <- terms(formula)
  event.history <- EHF$event.history
  n <- NROW(event.history)
  # }}}
  # {{{ time-constant design matrix Z
  ## Specials <- attr(Terms,"stripped.specials")
  ## strippedArguments <- attr(Terms,"stripped.arguments")
  Z <- EHF$prop
  if (!is.null(Z)) {
    power.args <- attr(EHF$prop, "arguments.terms")$power
    stopifnot(all(match(colnames(Z), names(power.args), nomatch = 0)))
    timeconst.power <- as.numeric(power.args[colnames(Z)])
    factorLevelsZ <- attr(EHF$prop, "levels")
    refLevelsZ <- lapply(factorLevelsZ, function(x) x[1])
    colnamesZ <- colnames(Z)
    dimZ <- NCOL(Z)
    fixed <- 1
    names(timeconst.power) <- colnames(Z)
    stopifnot(length(timeconst.power) == dimZ)
    if (!(all(timeconst.power %in% 0:2))) {
      stop("Only powers of time in 0,1,2 can be multipled to constant covariates")
    }
  } else {
    Z <- matrix(0, n, 1)
    dimZ <- 1
    colnamesZ <- NULL
    fixed <- 0
    factorLevelsZ <- NULL
    refLevelsZ <- NULL
    timeconst.power <- NULL
  }
  # }}}
  # {{{ time-varying design matrix X
  X <- EHF$timevar
  if (!is.null(X)) {
    test.args <- attr(EHF$timevar, "arguments.terms")$test
    stopifnot(all(match(colnames(X), names(test.args), nomatch = 0)))
    timevar.test <- as.numeric(test.args[colnames(X)])
    names(timevar.test) <- colnames(X)
    factorLevelsX <- attr(EHF$timevar, "levels")
    refLevelsX <- lapply(factorLevelsX, function(x) x[1])
    ## intercept
    X <- cbind("Intercept" = rep(1, n), X)
    timevar.test <- c(0, as.numeric(timevar.test))
    dimX <- NCOL(X)
  } else {
    ## X <- matrix(0,n,1)
    dimX <- 1
    colnamesX <- NULL
    factorLevelsX <- NULL
    refLevelsX <- NULL
    timevar.test <- 0
    ## intercept
    X <- cbind("Intercept" = rep(1, n))
  }
  stopifnot(length(timevar.test) == dimX)
  colnamesX <- colnames(X)
  if (!(all(timevar.test %in% 0:2))) {
    stop("Time power tests only available for powers 0,1,2")
  }
  theData <- data.frame(cbind(unclass(event.history), do.call("cbind", EHF[-1])))
  # }}}
  # {{{ event.history and order the data
  delayed <- !(is.null(attr(event.history, "entry.type"))) && !(attr(EHF$event.history, "entry.type") == "")
  if (delayed) {
    stop("Delayed entry is not (not yet) supported.")
  }
  cens.code <- attr(event.history, "cens.code")
  model.type <- attr(event.history, "model")
  states <- prodlim::getStates(event.history)
  stopifnot(model.type %in% c("survival", "competing.risks"))
  cens.type <- attr(event.history, "cens.type")
  stopifnot(cens.type %in% c("rightCensored", "uncensored"))
  neworder <- order(event.history[, "time"], -event.history[, "status"])
  event.history <- event.history[neworder, , drop = FALSE]
  Z <- Z[neworder, , drop = FALSE]
  X <- X[neworder, , drop = FALSE]
  theData <- theData[neworder, ]
  if (model.type[[1]] != "survival" && !("event" %in% colnames(event.history))) {
    warning("Only one cause of failure found in data.")
  }
  eventtime <- as.vector(event.history[, "time"])
  # time <- numeric(length(eventtime))
  if (model.type %in% c("competing.risks", "survival")) {
    if (cens.type %in% c("rightCensored", "uncensored")) {
      delta <- as.vector(event.history[, "status"])
      if (model.type == "competing.risks") {
        if (missing(cause)) {
          cause <- states[1]
          warning(paste("Argument cause is missing, analysing cause 1: ", cause, ". Other causes are:", paste(states[-1], collapse = ","), sep = ""))
        } else {
          if (!(cause %in% states)) stop(paste("Cause", cause, " is not among the causes in data; these are:", paste(states, collapse = ",")))
          cause <- match(cause, states, nomatch = 0)
        }
        ## event is 1 if the event of interest occured and 0 otherwise
        event <- event.history[, "event"] == cause
        if (sum(event) == 0) stop(paste("No events of type:", cause, "in data."))
      } else {
        event <- delta
      }
    } else {
      stop("Works only for right-censored data")
    }
  } else {
    stop("Response is neither competing risks nor survival.")
  }
  # }}}
  # {{{ cluster variable
  clusters <- EHF$cluster
  if (is.null(clusters)) {
    clusters <- 0:(NROW(X) - 1)
    # antclust <- NROW(X)
  } else {
    clusters <- as.integer(factor(clusters)) - 1
    # antclust <- length(unique(clusters))
  }
  # }}}
  # {{{ time points for timevarametric components
  if (missing(times) || length(times) == 0) {
    times <- NULL
  } else {
    times <- sort(unique(times))
  }
  # ntimes <- length(times)
  # }}}
  # {{{ ipcw model
  if (missing(cens.model)) cens.model <- "KM"
  if (length(grep("^km|^kaplan|^marg", cens.model, ignore.case = TRUE)) > 0) {
    cens.model <- "KM"
  } else {
    cens.model <- "cox"
  }
  # }}}
  enames <- colnames(event.history)
  if ("event" %in% enames) { # competing risks
    event.history[event.history[, "status"] == 0, "event"] <- 0
    event.history <- event.history[, -match("status", enames)]
  } else {
    colnames(event.history) <- sub("status", "event", colnames(event.history))
  }
  if ("entry" %in% enames) { # delayed entry
    timeregformula <- "timereg::Event(entry,time,event)"
  } else {
    timeregformula <- "timereg::Event(time,event)"
  }
  if ((ipos <- match("Intercept", colnamesX, nomatch = 0)) > 0) {
    timeregformula <- paste0(timeregformula, "~+1")
  } else {
    timeregformula <- paste0(timeregformula, "~-1")
  }
  if (length(colnamesZ) > 0) {
    # const <- timereg::const
    timeregformula <- paste0(
      timeregformula,
      "+",
      paste0(sapply(colnamesZ, function(z) {
        paste0("const(", z, ")")
      }), collapse = "+")
    )
  }
  if (length(colnamesX[-ipos]) > 0) {
    trdat <- data.frame(event.history, Z, X[, -ipos, drop = FALSE])
    timeregformula <- paste0(timeregformula, "+", paste0(colnamesX[-ipos], collapse = "+"))
  } else {
    trdat <- data.frame(event.history, Z)
  }
  # Surv <- survival::Surv
  out <- timereg::comp.risk(as.formula(timeregformula),
    data = trdat,
    model = trans,
    time.pow = timeconst.power,
    time.pow.test = timevar.test,
    times = times,
    cause = cause, ...
  )
  # }}}
  # {{{ prepare the output

  if (fixed == 1) {
    timeConstantCoef <- c(out$gamma)
    names(timeConstantCoef) <- colnamesZ
    timeConstantVar <- matrix(out$var.gamma, dimZ, dimZ, dimnames = list(colnamesZ, colnamesZ))
  } else {
    timeConstantCoef <- NULL
    timeConstantVar <- NULL
  }
  timeVaryingCoef <- out$cum
  timeVaryingVar <- out$var.cum
  score <- out$score
  timeConstantEffects <- list(coef = timeConstantCoef, var = timeConstantVar)
  class(timeConstantEffects) <- "timeConstantEffects"
  timeVaryingEffects <- list(
    coef = timeVaryingCoef,
    var = timeVaryingVar
  )
  attr(timeConstantEffects, "variables") <- names(factorLevelsZ)
  attr(timeVaryingEffects, "variables") <- names(factorLevelsX)
  class(timeVaryingEffects) <- "timeVaryingEffects"
  out <- list(
    call = call,
    response = event.history,
    design = list(
      Terms = Terms,
      const = colnamesZ,
      timevar = colnamesX,
      timepower = timeconst.power
    ),
    link = link,
    time = times,
    timeConstantEffects = timeConstantEffects,
    timeconst.power = timeconst.power,
    timeVaryingEffects = timeVaryingEffects,
    score = score,
    cens.model = cens.model,
    factorLevels = c(factorLevelsX, factorLevelsZ),
    refLevels = c(refLevelsX, refLevelsZ),
    "na.action" = attr(EHF, "na.action")
  )
  if (is.null(out$call$cause)) {
    out$call$cause <- cause
  }
  class(out) <- "riskRegression"
  return(out)
  # }}}
}
