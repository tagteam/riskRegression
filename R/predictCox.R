## * predictCox (documentation)
#' @title Fast computation of survival probabilities, hazards and cumulative hazards from Cox regression models
#' @name predictCox
#'
#' @description Fast routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a \code{survival::coxph} or \code{rms::cph} object
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata [data.frame or data.table]  Contain the values of the predictor variables
#' defining subject specific predictions.
#' Should have the same structure as the data set used to fit the \code{object}.
#' @param times [numeric vector] Time points at which to return
#' the estimated hazard/cumulative hazard/survival.
#' @param centered [logical] If \code{TRUE} return prediction at the
#'     mean values of the covariates \code{fit$mean}, if \code{FALSE}
#'     return a prediction for all covariates equal to zero.  in the
#'     linear predictor. Will be ignored if argument \code{newdata} is
#'     used. For internal use.
#' @param type [character vector] the type of predicted value. Choices are \itemize{
#'     \item \code{"hazard"} the baseline hazard function when
#'     argument \code{newdata} is not used and the hazard function
#'     when argument \code{newdata} is used.  \item \code{"cumhazard"}
#'     the cumulative baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  \item \code{"survival"}
#'     the survival baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  } Several choices can be
#'     combined in a vector of strings that match (no matter the case)
#'     strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param keep.strata [logical] If \code{TRUE} add the (newdata) strata
#'     to the output. Only if there any.
#' @param keep.times [logical] If \code{TRUE} add the evaluation times
#'     to the output.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the
#'     covariates used to make the prediction in the output list.
#' @param keep.infoVar [logical] For internal use.
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param diag [logical] when \code{FALSE} the hazard/cumlative hazard/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param average.iid [logical] If \code{TRUE} add the average of the influence function over \code{newdata} to the output.
#' @param store.iid [character] Implementation used to estimate the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}.
#'
#' @details
#' When the argument \code{newdata} is not specified, the function computes the baseline hazard estimate.
#' See (Ozenne et al., 2017) section "Handling of tied event times".
#'
#' Otherwise the function computes survival probabilities with confidence intervals/bands.
#' See (Ozenne et al., 2017) section "Confidence intervals and confidence bands for survival probabilities".
#' The survival is computed using the exponential approximation (equation 3).
#'
#' A detailed explanation about the meaning of the argument \code{store.iid} can be found
#' in (Ozenne et al., 2017) Appendix B "Saving the influence functions".
#'
#' The function is not compatible with time varying predictor variables.
#'
#' The centered argument enables us to reproduce the results obtained with the \code{basehaz}
#' function from the survival package but should not be modified by the user.
#'
#' The iid decomposition is output using an array containing the value of the influence
#' of each subject used to fit the object (dim 1),
#' for each subject in newdata (dim 3),
#' and each time (dim 2).
#'
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#'
#' @seealso
#' \code{\link{confint.predictCox}} to compute confidence intervals/bands.
#' \code{\link{autoplot.predictCox}} to display the predictions.

## * predictCox (examples)
#' @examples
##'
#' library(survival)
#' library(data.table)
#'
#' #### generate data ####
#' set.seed(10)
#' d <- sampleData(40, outcome = "survival") ## training dataset
#' nd <- sampleData(4, outcome = "survival") ## validation dataset
#' d$time <- round(d$time, 1) ## create tied events
#' # table(duplicated(d$time))
#'
#' #### stratified Cox model ####
#' fit <- coxph(Surv(time, event) ~ X1 + strata(X2) + X6,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' ## compute the baseline cumulative hazard
#' fit.haz <- predictCox(fit)
#' cbind(survival::basehaz(fit), fit.haz$cumhazard)
#'
#' ## compute individual specific cumulative hazard and survival probabilities
#' fit.pred <- predictCox(fit, newdata = nd, times = c(3, 8), se = TRUE, band = TRUE)
#' fit.pred
#'
#' ####  other examples ####
#' # one strata variable
#' fitS <- coxph(Surv(time, event) ~ strata(X1) + X2,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' predictCox(fitS)
#' predictCox(fitS, newdata = nd, times = 1)
#'
#' # two strata variables
#' set.seed(1)
#' d$U <- sample(letters[1:5], replace = TRUE, size = NROW(d))
#' d$V <- sample(letters[4:10], replace = TRUE, size = NROW(d))
#' nd$U <- sample(letters[1:5], replace = TRUE, size = NROW(nd))
#' nd$V <- sample(letters[4:10], replace = TRUE, size = NROW(nd))
#' fit2S <- coxph(Surv(time, event) ~ X1 + strata(U) + strata(V) + X2,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' cbind(survival::basehaz(fit2S), predictCox(fit2S, type = "cumhazard")$cumhazard)
#' predictCox(fit2S)
#' predictCox(fitS, newdata = nd, times = 3)
#'
#' # left truncation
#' test2 <- list(
#'   start = c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
#'   stop = c(2, 3, 6, 7, 8, 9, 9, 9, 14, 17),
#'   event = c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
#'   x = c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0)
#' )
#' m.cph <- coxph(Surv(start, stop, event) ~ 1, test2, x = TRUE)
#' as.data.table(predictCox(m.cph))
#'
#' basehaz(m.cph)
## * predictCox (code)
#' @rdname predictCox
#' @export
predictCox <- function(object,
                       times,
                       newdata = NULL,
                       centered = TRUE,
                       type = c("cumhazard", "survival"),
                       keep.strata = TRUE,
                       keep.times = TRUE,
                       keep.newdata = FALSE,
                       keep.infoVar = FALSE,
                       se = FALSE,
                       band = FALSE,
                       iid = FALSE,
                       confint = (se + band) > 0,
                       diag = FALSE,
                       average.iid = FALSE,
                       store.iid = "full") {
  call <- match.call()
  newdata <- copy(newdata)
  ## centering
  if (!is.null(newdata)) {
    if (inherits(centered, "data.frame")) {
      df.reference <- centered
      centered2 <- TRUE ## for the linear predictor of the hazard
    } else {
      df.reference <- NULL
      centered2 <- centered ## for the linear predictor of the hazard
    }
    centered <- TRUE ## for the linear predictor of the baseline hazard
  } else {
    centered2 <- FALSE
  }

  ## ** Extract elements from object
  if (missing(times)) {
    nTimes <- 0
    times <- numeric(0)
  } else {
    nTimes <- length(times)
  }
  needOrder <- (nTimes[1] > 0 && is.unsorted(times))
  if (all(!is.na(times)) && needOrder) {
    order.times <- order(times)
    oorder.times <- order(order.times)
    times.sorted <- sort(times)
  } else {
    if (nTimes == 0) {
      times.sorted <- numeric(0)
    } else {
      times.sorted <- times
      order.times <- 1:nTimes
      oorder.times <- 1:nTimes
    }
  }
  object.n <- coxN(object)
  object.modelFrame <- coxModelFrame(object)
  infoVar <- coxVariableName(object, model.frame = object.modelFrame)
  object.baseEstimator <- coxBaseEstimator(object)

  ## ease access
  is.strata <- infoVar$is.strata
  object.levelStrata <- levels(object.modelFrame$strata) ## levels of the strata variable
  nStrata <- length(object.levelStrata) ## number of strata
  nVar.lp <- length(infoVar$lpvars) ## number of variables in the linear predictor

  ## ** normalize model frame
  ## convert strata to numeric
  object.modelFrame[, c("strata.num") := as.numeric(.SD$strata) - 1]

  ## linear predictor
  ## if we predict the hazard for newdata then there is no need to center the covariates
  object.modelFrame[, c("eXb") := exp(coxLP(object, data = NULL, center = if (is.null(newdata)) {
    centered
  } else {
    FALSE
  }))]
  ## add linear predictor and remove useless columns
  rm.name <- setdiff(names(object.modelFrame), c("start", "stop", "status", "eXb", "strata", "strata.num"))
  if (length(rm.name) > 0) {
    object.modelFrame[, c(rm.name) := NULL]
  }

  ## sort the data
  object.modelFrame[, c("statusM1") := 1 - .SD$status] ## sort by statusM1 such that deaths appear first and then censored events
  object.modelFrame[, c("XXXindexXXX") := 1:.N] ## keep track of the initial positions (useful when calling calcSeCox)
  data.table::setkeyv(object.modelFrame, c("strata.num", "stop", "start", "statusM1"))

  ## last event time in each strata
  if (!is.null(attr(times, "etimes.max"))) { ## specified by the user
    etimes.max <- attr(times, "etimes.max")
    attr(times, "etimes.max") <- NULL
    attr(times.sorted, "etimes.max") <- etimes.max
  } else if (is.strata) { ## based on the data
    if (nVar.lp == 0) {
      iDTtempo <- object.modelFrame[, .SD[which.max(.SD$stop)], by = "strata.num"]
      etimes.max <- iDTtempo[, if (.SD$status == 1) {
        1e12
      } else {
        .SD$stop
      }, by = "strata.num"][[2]]
    } else {
      etimes.max <- object.modelFrame[, max(.SD$stop), by = "strata.num"][[2]]
    }
  } else {
    if (nVar.lp == 0 && (utils::tail(object.modelFrame$status, 1) == 1)) { ## no covariates and ends by a death
      etimes.max <- 1e12
    } else {
      etimes.max <- max(object.modelFrame[["stop"]])
    }
  }

  ## ** checks
  ## check user imputs
  if (nTimes[1] > 0 && any(is.na(times))) {
    stop("Missing (NA) values in argument \'times\' are not allowed.\n")
  }
  type <- tolower(type)
  if (any(type %in% c("lp", "hazard", "cumhazard", "survival") == FALSE)) {
    stop("type can only be \"lp\", \"hazard\", \"cumhazard\" or/and \"survival\" \n")
  }
  if (is.null(newdata) && "lp" %in% type) {
    stop("Cannot evaluate the linear predictor when argument \'newdata\' is missing. \n")
  }
  if (length(times) > 1 && "lp" %in% type) {
    stop("Cannot evaluate the linear predictor when there are multiple timepoints. \n")
  }
  ## predictCox is not compatible with all coxph/cph object (i.e. only handle only simple cox models)
  if (!is.null(object$weights) && !all(object$weights == 1)) {
    stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if (!is.null(object$naive.var)) {
    stop("predictCox does not know how to handle frailty.")
  }
  ## if(any(object.modelFrame[["start"]]!=0)){
  ##     warning("The current version of predictCox was not designed to handle left censoring \n",
  ##             "The function may be used on own risks \n")
  ## }
  if (object.baseEstimator == "exact") {
    stop("Prediction with exact handling of ties is not implemented.\n")
  }
  if (!is.null(object$call$tt)) {
    stop("predictCox does not know how to handle time varying effects.\n")
  }
  ## convergence issue
  if (!is.null(coef(object)) && any(is.na(coef(object)))) {
    print(coef(object))
    stop("One or several parameters of the regression model have no value, i.e., a value 'NA'.\n")
  }
  ## prediction
  if (missing(newdata) && (se || iid || average.iid)) {
    stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
  }
  if (!is.null(newdata)) {
    if (nTimes[1] == 0 && !identical(as.character(type), "lp")) {
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    if (!is.vector(times)) {
      stop("Argument \'times\' must be a vector \n")
    }
    name.regressor <- c(infoVar$lpvars.original, infoVar$stratavars.original)
    if (length(name.regressor) > 0 && any(name.regressor %in% names(newdata) == FALSE)) {
      stop(
        "Missing variables in argument \'newdata\': \"",
        paste0(setdiff(name.regressor, names(newdata)), collapse = "\" \""),
        "\"\n"
      )
    }
    if (se[1] && ("hazard" %in% type)) {
      stop("Standard error cannot be computed for the hazard \n")
    }
    if (band[1] && ("hazard" %in% type)) {
      stop("confidence bands cannot be computed for the hazard \n")
    }
  }
  ## diag argument
  if (!is.logical(diag)) {
    stop("Argument \'diag\' must be of type logical \n")
  }
  if (diag) {
    if (NROW(newdata) != length(times)) {
      stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
    }
  }
  if (average.iid == TRUE && !is.null(attr(average.iid, "factor"))) {
    if (is.null(store.iid) && !is.null(object$iid$store.iid)) {
      store.iid <- object$iid$store.iid
    }
    if (store.iid == "full") {
      if (iid) {
        stop(
          "Attribute \"factor\" of argument \'average.iid\' not available when \'iid\' is TRUE with argument \'store.iid\' set to \"full\" \n",
          "Consider setting  \'store.iid\' set to \"minimal\" \n"
        )
      }
      if (se) {
        stop(
          "Attribute \"factor\" of argument \'average.iid\' not available when \'se\' is TRUE with argument \'store.iid\' set to \"full\" \n",
          "Consider setting  \'store.iid\' set to \"minimal\" \n"
        )
      }
    }
    test.list <- !is.list(attr(average.iid, "factor"))
    if (test.list) {
      stop("Attribute \"factor\" of argument \'average.iid\' must be a list \n")
    }
    test.matrix <- any(unlist(lapply(attr(average.iid, "factor"), is.matrix)) == FALSE)
    if (test.matrix) {
      stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices \n")
    }
    for (iFactor in 1:length(attr(average.iid, "factor"))) { ## iFactor <- 1
      ## when only one column and diag = FALSE, use the same weights at all times
      if ((diag == FALSE) && (NCOL(attr(average.iid, "factor")[[iFactor]]) == 1) && (nTimes > 1)) {
        attr(average.iid, "factor")[[iFactor]] <- matrix(attr(average.iid, "factor")[[iFactor]][, 1],
          nrow = NROW(attr(average.iid, "factor")[[iFactor]]),
          ncol = nTimes, byrow = FALSE
        )
      }
      ## check dimensions
      if (any(dim(attr(average.iid, "factor")[[iFactor]]) != c(NROW(newdata), diag + (1 - diag) * nTimes))) {
        stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ", new.n, ",", diag + (1 - diag) * nTimes, " \n")
      }
    }
  }

  ## ** baseline hazard

  ## compute the baseline hazard
  Lambda0 <- baseHaz_cpp(
    starttimes = object.modelFrame$start,
    stoptimes = object.modelFrame$stop,
    status = object.modelFrame$status,
    eXb = object.modelFrame$eXb,
    strata = object.modelFrame$strata.num,
    nPatients = object.n,
    nStrata = nStrata,
    emaxtimes = etimes.max,
    predtimes = times.sorted,
    cause = 1,
    Efron = (object.baseEstimator == "efron")
  )
  ## restaure strata levels
  Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata - 1), labels = object.levelStrata)

  ## ** compute cumlative hazard and survival
  if (is.null(newdata)) {
    if (!("hazard" %in% type)) {
      Lambda0$hazard <- NULL
    }
    if ("survival" %in% type) { ## must be before cumhazard
      Lambda0$survival <- exp(-Lambda0$cumhazard)
    }
    if (!("cumhazard" %in% type)) {
      Lambda0$cumhazard <- NULL
    }
    if (keep.times == FALSE) {
      Lambda0$times <- NULL
    }
    if (keep.strata[[1]] == FALSE || (is.null(call$keep.strata) && !is.strata)) {
      Lambda0$strata <- NULL
    }
    if (keep.newdata[1] == TRUE) {
      Lambda0$newdata <- object.modelFrame
    }
    add.list <- list(
      lastEventTime = etimes.max,
      se = FALSE,
      band = FALSE,
      type = type,
      nTimes = nTimes,
      baseline = TRUE,
      var.lp = infoVar$lpvars.original,
      var.strata = infoVar$stratavars.original
    )
    if (keep.infoVar) {
      add.list$infoVar <- infoVar
    }

    Lambda0[names(add.list)] <- add.list
    class(Lambda0) <- "predictCox"
    return(Lambda0)
  } else {
    out <- list()
    ## *** reformat newdata (compute linear predictor and strata)
    new.n <- NROW(newdata)
    setDT(newdata)

    Xb <- coxLP(object, data = newdata, center = FALSE)
    if ("lp" %in% type) {
      out$lp <- cbind(Xb)
      lp.iid <- centered2 ## when ask for centered then we need the iid for exporting the correct se
    } else {
      lp.iid <- FALSE
    }
    new.eXb <- exp(Xb)
    new.strata <- coxStrata(object,
      data = newdata,
      sterms = infoVar$strata.sterms,
      strata.vars = infoVar$stratavars,
      strata.levels = infoVar$strata.levels
    )

    new.levelStrata <- levels(droplevels(new.strata))

    ## *** subject specific hazard
    if (is.strata == FALSE && !identical(as.character(type), "lp")) {
      if (diag) {
        if (needOrder) {
          iTimes <- prodlim::sindex(jump.times = Lambda0$times, eval.times = times.sorted[oorder.times])
        } else {
          iTimes <- prodlim::sindex(jump.times = Lambda0$times, eval.times = times.sorted)
        }
      }

      if ("hazard" %in% type) {
        if (diag) {
          out$hazard <- cbind(new.eXb * Lambda0$hazard[iTimes])
        } else {
          out$hazard <- (new.eXb %o% Lambda0$hazard)
          if (needOrder) out$hazard <- out$hazard[, oorder.times, drop = 0L]
        }
      }
      if ("cumhazard" %in% type || "survival" %in% type) {
        if (diag) {
          cumhazard <- cbind(new.eXb * Lambda0$cumhazard[iTimes])
        } else {
          cumhazard <- new.eXb %o% Lambda0$cumhazard
          if (needOrder) {
            cumhazard <- cumhazard[, oorder.times, drop = 0L]
          }
        }
        if ("cumhazard" %in% type) {
          out$cumhazard <- cumhazard
        }
        if ("survival" %in% type) {
          out$survival <- exp(-cumhazard)
        }
      }
    } else if (!identical(as.character(type), "lp")) {
      ## initialization
      if ("hazard" %in% type) {
        out$hazard <- matrix(0, nrow = new.n, ncol = nTimes * (1 - diag) + diag)
      }
      if ("cumhazard" %in% type) {
        out$cumhazard <- matrix(NA, nrow = new.n, ncol = nTimes * (1 - diag) + diag)
      }
      if ("survival" %in% type) {
        out$survival <- matrix(NA, nrow = new.n, ncol = nTimes * (1 - diag) + diag)
      }

      ## loop across strata
      for (S in new.levelStrata) { ## S <- 1
        id.S <- which(Lambda0$strata == S)
        newid.S <- which(new.strata == S)

        if (diag) {
          if (needOrder) {
            iSTimes <- prodlim::sindex(jump.times = Lambda0$times[id.S], eval.times = times.sorted[oorder.times[newid.S]])
          } else {
            iSTimes <- prodlim::sindex(jump.times = Lambda0$times[id.S], eval.times = times.sorted[newid.S])
          }
        }

        if ("hazard" %in% type) {
          if (diag) {
            out$hazard[newid.S] <- new.eXb[newid.S] * Lambda0$hazard[id.S][iSTimes]
          } else {
            out$hazard[newid.S, ] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S]
            if (needOrder) {
              out$hazard[newid.S, ] <- out$hazard[newid.S, oorder.times, drop = 0L]
            }
          }
        }
        if ("cumhazard" %in% type || "survival" %in% type) {
          if (diag) {
            cumhazard.S <- cbind(new.eXb[newid.S] * Lambda0$cumhazard[id.S][iSTimes])
          } else {
            cumhazard.S <- new.eXb[newid.S] %o% Lambda0$cumhazard[id.S]
            if (needOrder) {
              cumhazard.S <- cumhazard.S[, oorder.times, drop = 0L]
            }
          }

          if ("cumhazard" %in% type) {
            out$cumhazard[newid.S, ] <- cumhazard.S
          }
          if ("survival" %in% type) {
            out$survival[newid.S, ] <- exp(-cumhazard.S)
          }
        }
      }
    }

    if (se[[1]] || band[[1]] || iid[[1]] || average.iid[[1]]) {
      if (nVar.lp > 0) {
        ## get the (new) design matrix
        new.LPdata <- model.matrix(object, data = newdata)
        if (NROW(new.LPdata) != NROW(newdata)) {
          stop(
            "NROW of the design matrix and newdata differ. \n",
            "Maybe because newdata contains NA values \n"
          )
        }
        if (any(sort(colnames(new.LPdata)) != sort(names(coef(object))))) {
          stop(
            "Names of the design matrix and model parameters differ. \n",
            "Possible error in model.matrix due to special operator in the formula. \n"
          )
        }
      } else {
        new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
      }

      ## restaure original ordering
      data.table::setkeyv(object.modelFrame, "XXXindexXXX")
      if (diag) {
        Lambda0$oorder.times <- oorder.times
      } else {
        Lambda0$oorder.times <- 1:nTimes
      }

      ## Computation of the influence function and/or the standard error
      export <- c("iid"[(iid + band + lp.iid) > 0], "se"[(se + band) > 0], "average.iid"[average.iid == TRUE])
      if (!is.null(attr(average.iid, "factor"))) {
        if (diag) {
          attr(export, "factor") <- attr(average.iid, "factor")
        } else {
          ## re-order columns according to times
          attr(export, "factor") <- lapply(attr(average.iid, "factor"), function(iF) {
            iF[, order.times, drop = FALSE]
          })
        }
      }
      if (diag) {
        times2 <- times
      } else {
        times2 <- times.sorted
      }
      attr(times2, "etimes.max") <- attr(times.sorted, "etimes.max")

      outSE <- calcSeCox(object,
        times = times2,
        nTimes = nTimes,
        type = type,
        diag = diag,
        Lambda0 = Lambda0,
        object.n = object.n,
        object.time = object.modelFrame$stop,
        object.eXb = object.modelFrame$eXb,
        object.strata = object.modelFrame$strata,
        nStrata = nStrata,
        new.n = new.n,
        new.eXb = new.eXb,
        new.LPdata = new.LPdata,
        new.strata = new.strata,
        new.survival = if (diag) {
          out$survival
        } else {
          out$survival[, order.times, drop = FALSE]
        },
        nVar.lp = nVar.lp,
        export = export,
        store.iid = store.iid
      )

      ## restaure orginal time ordering
      if ((iid + band) > 0) {
        if ("lp" %in% type) {
          out$lp.iid <- outSE$lp.iid
        }
        if ("hazard" %in% type) {
          if (needOrder[1] && (diag[1] == FALSE)) {
            out$hazard.iid <- outSE$hazard.iid[, oorder.times, , drop = 0L]
          } else {
            out$hazard.iid <- outSE$hazard.iid
          }
        }
        if ("cumhazard" %in% type) {
          if (needOrder[1] && (diag[1] == FALSE)) {
            out$cumhazard.iid <- outSE$cumhazard.iid[, oorder.times, , drop = 0L]
          } else {
            out$cumhazard.iid <- outSE$cumhazard.iid
          }
        }
        if ("survival" %in% type) {
          if (needOrder[1] && (diag[1] == FALSE)) {
            out$survival.iid <- outSE$survival.iid[, oorder.times, , drop = 0L]
          } else {
            out$survival.iid <- outSE$survival.iid
          }
        }
      }
      if (average.iid == TRUE) {
        if ("lp" %in% type) {
          out$lp.average.iid <- outSE$lp.average.iid
        }
        if ("hazard" %in% type) {
          if (needOrder && (diag[1] == FALSE)) {
            if (is.list(outSE$hazard.average.iid)) {
              out$hazard.average.iid <- lapply(outSE$hazard.average.iid, function(iIID) {
                iIID[, oorder.times, drop = 0L]
              })
            } else {
              out$hazard.average.iid <- outSE$hazard.average.iid[, oorder.times, drop = 0L]
            }
          } else {
            out$hazard.average.iid <- outSE$hazard.average.iid
          }
        }
        if ("cumhazard" %in% type) {
          if (needOrder && (diag[1] == FALSE)) {
            if (is.list(outSE$cumhazard.average.iid)) {
              out$cumhazard.average.iid <- lapply(outSE$cumhazard.average.iid, function(iIID) {
                iIID[, oorder.times, drop = 0L]
              })
            } else {
              out$cumhazard.average.iid <- outSE$cumhazard.average.iid[, oorder.times, drop = 0L]
            }
          } else {
            out$cumhazard.average.iid <- outSE$cumhazard.average.iid
          }
        }
        if ("survival" %in% type) {
          if (needOrder && (diag[1] == FALSE)) {
            if (is.list(outSE$survival.average.iid)) {
              out$survival.average.iid <- lapply(outSE$survival.average.iid, function(iIID) {
                iIID[, oorder.times, drop = 0L]
              })
            } else {
              out$survival.average.iid <- outSE$survival.average.iid[, oorder.times, drop = 0L]
            }
          } else {
            out$survival.average.iid <- outSE$survival.average.iid
          }
        }

        if (is.list(attr(average.iid, "factor"))) {
          if ("hazard" %in% type) {
            names(out$hazard.average.iid) <- names(attr(average.iid, "factor"))
          }
          if ("cumhazard" %in% type) {
            names(out$cumhazard.average.iid) <- names(attr(average.iid, "factor"))
          }
          if ("survival" %in% type) {
            names(out$survival.average.iid) <- names(attr(average.iid, "factor"))
          }
        }
      }
      if ((se + band) > 0) {
        if ("lp" %in% type) {
          out$lp.se <- outSE$lp.se
        }
        if ("cumhazard" %in% type) {
          if (needOrder && (diag[1] == FALSE)) {
            out$cumhazard.se <- outSE$cumhazard.se[, oorder.times, drop = 0L]
          } else {
            out$cumhazard.se <- outSE$cumhazard.se
          }
        }
        if ("survival" %in% type) {
          if (needOrder && (diag[1] == FALSE)) {
            out$survival.se <- outSE$survival.se[, oorder.times, drop = 0L]
          } else {
            out$survival.se <- outSE$survival.se
          }
        }
      }
    }

    ## ** substract reference
    if ("lp" %in% type && centered2) {
      if (is.null(df.reference)) {
        data <- try(eval(object$call$data), silent = TRUE)
        if (inherits(x = data, what = "try-error")) {
          stop(
            "Could not evaluate the dataset used to fit the model to define a reference level. \n",
            "Set argument \'centered\' to FALSE or to a data.frame definining the reference level. \n"
          )
        }
        var.original <- infoVar$lpvars.original

        ls.ref <- lapply(var.original, function(iVar) {
          if (is.numeric(data[[iVar]])) {
            return(unname(mean(data[[iVar]])))
          } else if (is.factor(data[[iVar]])) {
            return(unname(factor(levels(data[[iVar]])[1], levels(data[[iVar]]))))
          } else if (is.character(data[[iVar]])) {
            return(unname(sort(unique(data[[iVar]])[1])))
          }
        })
        df.reference <- as.data.frame(setNames(ls.ref, var.original))
      }

      ls.args <- as.list(call)[-1]
      ls.args$newdata <- df.reference
      ls.args$centered <- FALSE
      ls.args$type <- "lp"
      ls.args$se <- (se + band > 0)
      ls.args$iid <- (iid + se + band > 0)
      ls.args$band <- FALSE
      ls.args$confint <- FALSE
      outRef <- do.call(predictCox, args = ls.args)

      out$lp <- out$lp - as.double(outRef$lp)
      if (band[1] || se[1]) {
        out$lp.se <- cbind(sqrt(colSums(colCenter_cpp(outSE$lp.iid, outRef$lp.iid)^2))) ## use outSe$lp.iid instead of out$lp.iid for when not exporting the iid
      }
      if (band[1] || iid[1]) {
        out$lp.iid <- colCenter_cpp(out$lp.iid, outRef$lp.iid)
      }
    }

    ## ** add information to the predictions
    add.list <- list(
      lastEventTime = etimes.max,
      se = se,
      band = band,
      type = type,
      diag = diag,
      nTimes = nTimes,
      baseline = FALSE,
      var.lp = infoVar$lpvars.original,
      var.strata = infoVar$stratavars.original
    )
    if (keep.times == TRUE) {
      add.list$times <- times
    }
    if (is.strata[1] && keep.strata[1] == TRUE) {
      add.list$strata <- new.strata
    }

    if (keep.infoVar) {
      add.list$infoVar <- infoVar
    }
    all.covars <- c(infoVar$stratavars.original, infoVar$lpvars.original)
    if (keep.newdata[1] == TRUE && length(all.covars) > 0) {
      if (data.table::is.data.table(newdata)) {
        add.list$newdata <- newdata[, all.covars, with = FALSE]
      } else {
        add.list$newdata <- newdata[, all.covars, drop = FALSE]
      }
    }
    out[names(add.list)] <- add.list
    class(out) <- "predictCox"

    ## ** confidence intervals/bands
    if (confint) {
      out <- stats::confint(out)
    }
    if (band[1] && se[1] == FALSE) {
      out[paste0(type, ".se")] <- NULL
    }
    if (band[1] && iid[1] == FALSE) {
      out[paste0(type, ".iid")] <- NULL
    }

    return(out)
  }
}
