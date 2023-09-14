#' Cause-specific Cox proportional hazard regression
#'
#' Interface for fitting cause-specific Cox proportional hazard regression
#' models in competing risk.
#'
#' The causes and their order
#' are determined by \code{prodlim::getStates()} applied to the Hist object.
#' @param formula Either a single \code{Hist} formula or a list of formulas.
#' If it is a list it must contain as many \code{Hist} formulas as there are
#' causes when \code{surv.type="hazard"} and exactly two formulas when \code{surv.type="survival"}.
#' If it is a list the first formula is used for the cause of interest specific Cox regression
#' and the other formula(s) either for the other cause specific Cox regression(s) or for the
#' Cox regression of the combined event where each cause counts as event. Note that when only one
#' formula is given the covariates enter in exactly the same way into all Cox regression analyses.
#' @param data A data in which to fit the models.
#' @param cause The cause of interest. Defaults to the first cause (see Details).
#' @param surv.type Either \code{"hazard"} (the default) or
#' \code{"survival"}.  If \code{"hazard"} fit cause-specific Cox
#' regression models for all causes.  If \code{"survival"} fit one
#' cause-specific Cox regression model for the cause of interest and
#' also a Cox regression model for event-free survival.
#' @param fitter Routine to fit the Cox regression models.
#' If \code{coxph} use \code{survival::coxph} else use \code{rms::cph}.
#' @param ... Arguments given to \code{fitter}, e.g., \code{coxph}.
#' @return \item{models }{a list with the fitted (cause-specific) Cox
#' regression objects} \item{response }{the event history response }
#' \item{eventTimes }{the sorted (unique) event times } \item{surv.type }{the
#' value of \code{surv.type}} \item{theCause }{the cause of interest. see
#' \code{cause}} \item{causes }{the other causes} %% ...
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk} and Ulla B. Mogensen
#' @references
#' B. Ozenne, A. L. Soerensen, T.H. Scheike, C.T. Torp-Pedersen,
#' and T.A. Gerds. riskregression: Predicting the risk
#' of an event using Cox regression models. R Journal, 9(2):440--460, 2017.
#'
#' J Benichou and Mitchell H Gail. Estimates of absolute cause-specific risk
#' in cohort studies. Biometrics, pages 813--826, 1990.
#'
#' T.A. Gerds, T.H. Scheike, and P.K. Andersen. Absolute risk regression for
#' competing risks: Interpretation, link functions, and prediction. Statistics
#' in Medicine, 31(29):3921--3930, 2012.
#'
#' @seealso \code{\link{coxph}}
#' @keywords survival
#' @examples
#'
#' #' library(prodlim)
#' #' library(survival)
#' #' data(Melanoma)
#' #' ## fit two cause-specific Cox models
#' #' ## different formula for the two causes
#' #' fit1 <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
#' #'             data=Melanoma)
#' #' print(fit1)
#' #'
#' \dontrun{
#' #' library(Publish)
#' #' publish(fit1)
#' #'
#' }
#' #'
#' #' ## model hazard of all cause mortality instead of hazard of type 2
#' #' fit1a <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
#' #'              data=Melanoma,
#' #'              surv.type="surv")
#' #'
#' #' ## the predicted probabilities are similar
#' #' plot(predictRisk(fit1,times=500,cause=1,newdata=Melanoma),
#' #'      predictRisk(fit1a,times=500,cause=1,newdata=Melanoma))
#' #'
#' #' ## special case where cause 2 has no covariates
#' #' fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~1),
#' #'              data=Melanoma)
#' #' print(fit1b)
#' #' predict(fit1b,cause=1,times=100,newdata=Melanoma)
#' #'
#' #'
#' #' ## same formula for both causes
#' #' fit2 <- CSC(Hist(time,status)~invasion+epicel+age,
#' #'             data=Melanoma)
#' #' print(fit2)
#' #'
#' #' ## combine a cause-specific Cox regression model for cause 2
#' #' ## and a Cox regression model for the event-free survival:
#' #' ## different formula for cause 2 and event-free survival
#' #' fit3 <- CSC(list(Hist(time,status)~sex+invasion+epicel+age,
#' #'                  Hist(time,status)~invasion+epicel+age),
#' #'             surv.type="surv",
#' #'             data=Melanoma)
#' #' print(fit3)
#' #'
#' #' ## same formula for both causes
#' #' fit4 <- CSC(Hist(time,status)~invasion+epicel+age,
#' #'             data=Melanoma,
#' #'             surv.type="surv")
#' #' print(fit4)
#' #'
#' #' ## strata
#' #' fit5 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
#' #'             data=Melanoma,
#' #'             surv.type="surv")
#' #' print(fit5)
#' #'
#' #' ## sanity checks
#' #'
#' #' cox1 <- coxph(Surv(time,status==1)~invasion+epicel+age+strata(sex),data=Melanoma)
#' #' cox2 <- coxph(Surv(time,status!=0)~invasion+epicel+age+strata(sex),data=Melanoma)
#' #' all.equal(coef(cox1),coef(fit5$models[[1]]))
#' #' all.equal(coef(cox2),coef(fit5$models[[2]]))
#' #'
#' #' ## predictions
#' #' ##
#' #' ## surv.type = "hazard": predictions for both causes can be extracted
#' #' ## from the same fit
#' #' fit2 <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma)
#' #' predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
#' #' predictRisk(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
#' #' predictRisk(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
#' #' predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
#' #' predict(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
#' #'
#' #' ## surv.type = "surv" we need to change the cause of interest
#' #' library(survival)
#' #' fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
#' #'             data=Melanoma,
#' #'             surv.type="surv",cause=2)
#' #' ## now this does not work
#' #' try(predictRisk(fit5.2,cause=1,newdata=Melanoma,times=4))
#' #'
#' #' ## but this does
#' #' predictRisk(fit5.2,cause=2,newdata=Melanoma,times=100)
#' #' predict(fit5.2,cause=2,newdata=Melanoma,times=100)
#' #' predict(fit5.2,cause=2,newdata=Melanoma[4,],times=100)
#' #'
#' @export
CSC <- function(formula,
                data,
                cause,
                surv.type = "hazard",
                fitter = "coxph",
                ## strip.environment
                ...) {
  fitter <- match.arg(fitter, c("coxph", "cph", "phreg"))
  # {{{ type
  surv.type <- match.arg(surv.type, c("hazard", "survival"))
  # }}}
  # {{{ formulae & response
  if (inherits(x = formula, what = "formula")) formula <- list(formula)
  call <- match.call()
  # get outcome information from formula
  Rform <- update(formula[[1]], ".~1")
  ## mf <- stats::model.frame(Rform, data = data, na.action = na.omit)
  ## response <- stats::model.response(mf)
  response <- eval(Rform[[2]], envir = data)
  if (any(is.na(response))) {
    stop("Event history response may not contain missing values")
  }
  time <- response[, "time"]
  status <- response[, "status"]
  event <- prodlim::getEvent(response)
  if ("entry" %in% colnames(response)) {
    entry <- response[, "entry"]
  } else {
    if (fitter == "phreg") {
      entry <- 0
    } else {
      entry <- NULL
    }
  }
  if (any(entry > time)) stop("entry > time detected. Entry time into the study must be strictly greater than outcome time.")
  ## remove event history variables from data
  if (any((this <- match(all.vars(Rform), names(data), nomatch = 0)) > 0)) {
    if (data.table::is.data.table(data)) {
      data <- data[, -this, with = FALSE]
    } else {
      data <- data[, -this]
    }
  }
  # }}}
  # {{{ sorted unique event times
  eventTimes <- unique(sort(as.numeric(time[status != 0])))
  # }}}
  # {{{ causes
  causes <- prodlim::getStates(response)
  if (surv.type == "hazard") {
    NC <- length(causes)
  } else {
    NC <- 2
  } # cause of interest and overall survival
  if (length(formula) != NC[1] && length(formula) > 1) stop("Wrong number of formulae. Should be one for each cause ", NC, ".")
  if (length(formula) == 1) {
    ## warning("The same formula used for all causes")
    formula <- lapply(1:NC, function(x) formula[[1]])
  }
  # }}}
  # {{{ find the cause of interest
  if (missing(cause)) {
    theCause <- causes[1]
  } else {
    if ((foundCause <- match(as.character(cause), causes, nomatch = 0)) == 0) {
      stop(paste0(
        "Cannot find all requested cause(s) ...\n\n",
        "Requested cause(s): ", paste0(cause, collapse = ", "),
        "\n Available causes: ", paste(causes, collapse = ", "),
        "\n"
      ))
    } else {
      theCause <- causes[foundCause]
    }
  }
  otherCauses <- causes[-match(theCause, causes)]
  # }}}
  # {{{ fit Cox models
  if (surv.type == "hazard") {
    nmodels <- NC
  } else {
    nmodels <- 2
  }
  CoxModels <- lapply(1:nmodels, function(x) {
    if (surv.type == "hazard") {
      if (x == 1) {
        causeX <- theCause
      } else {
        causeX <- otherCauses[x - 1]
      }
    } else {
      causeX <- theCause
    }
    if (surv.type == "hazard") {
      statusX <- as.numeric(event == causeX)
    } else {
      if (x == 1) {
        statusX <- 1 * (event == causeX)
      } else { ## event-free status
        statusX <- response[, "status"]
      }
    }
    if (is.null(entry)) {
      workData <- data.frame(time = time, status = statusX)
    } else {
      workData <- data.frame(time = time, status = statusX, entry = entry)
    }
    if (any(this <- match(names(data), names(workData), nomatch = 0) > 0)) {
      warning(paste("Variables named", paste(names(data)[this], collapse = ", "), "in data will be ignored."))
      if (is.data.table(data)) {
        data <- data[, -this, with = FALSE]
      } else {
        data <- data[, -this, drop = FALSE]
      }
    }
    workData <- cbind(workData, data)
    if (is.null(entry)) {
      survresponse <- "survival::Surv(time, status)"
    } else {
      survresponse <- "survival::Surv(entry, time, status)"
    }
    ## check whether right hand side of formula includes ~.
    allvars <- all.vars(formula[[x]])
    if (any(grepl("^\\.$", allvars))) {
      formulaXX <- as.formula(paste0(survresponse, "~."))
    } else {
      formulaXX <- update(formula[[x]], paste0(survresponse, "~."))
    }
    # previous
    # formulaXX <- update(formula[[x]],paste0(survresponse,"~."))
    ## as.formula(paste(survresponse,
    ## as.character(delete.response(terms.formula(formulaX)))[[2]],
    ## sep="~"))
    args <- list(formulaXX, data = workData)
    extra.args <- list(...)
    if (fitter == "coxph") {
      fit <- do.call("coxph", c(args, list(x = TRUE, y = TRUE), extra.args))
    } else if (fitter == "cph") {
      fit <- do.call("cph", c(args, list(surv = TRUE, x = TRUE, y = TRUE), extra.args))
    } else if (fitter == "phreg") {
      fit <- do.call("phreg", c(args, extra.args))
    }
    ## fit$formula <- terms(fit$formula)
    ## fit$call$formula <- terms(formulaXX)
    ## fit$call$formula <- fit$formula
    ## fit$call$data <- NULL
    fit$call$formula <- formulaXX
    fit$call$data <- workData
    fit
  })
  if (surv.type == "hazard") {
    names(CoxModels) <- paste("Cause", c(theCause, otherCauses))
  } else {
    names(CoxModels) <- c(paste("Cause", theCause), "OverallSurvival")
  }
  # }}}
  out <- list(
    call = call,
    models = CoxModels,
    response = response,
    eventTimes = eventTimes,
    surv.type = surv.type,
    fitter = fitter,
    theCause = theCause,
    causes = c(theCause, otherCauses)
  )
  class(out) <- "CauseSpecificCox"
  out
}
