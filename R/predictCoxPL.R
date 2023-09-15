### predictCoxPL.R ---
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep  4 2017 (16:43)
## Version:
## last-updated: Dec 20 2021 (12:24)
##           By: Brice Ozenne
##     Update #: 172
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

## * predictCoxPL (documentation)
#' @title Computation of survival probabilities from Cox regression models using the product limit estimator.
#' @description Same as predictCox except that the survival is estimated using the product limit estimator.
#' @name predictCoxPL
#'
#' @inheritParams predictCox
#' @param ... additional arguments to be passed to \code{\link{predictCox}}.
#'
#' @details Note: the iid and standard errors are computed using the exponential approximation.

## * predictCoxPL (code)
#' @examples
#' library(survival)
#'
#' #### generate data ####
#' set.seed(10)
#' d <- sampleData(40, outcome = "survival")
#' nd <- sampleData(4, outcome = "survival")
#' d$time <- round(d$time, 1)
#'
#' #### Cox model ####
#' fit <- coxph(Surv(time, event) ~ X1 + X2 + X6,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' ## exponential approximation
#' predictCox(fit, newdata = d, times = 1:5)
#'
#' ## product limit
#' predictCoxPL(fit, newdata = d, times = 1:5)
#'
#' #### stratified Cox model ####
#' fitS <- coxph(Surv(time, event) ~ X1 + strata(X2) + X6,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' ## exponential approximation
#' predictCox(fitS, newdata = d, times = 1:5)
#'
#' ## product limit
#' predictCoxPL(fitS, newdata = d, times = 1:5)
#'
#' #### fully stratified Cox model ####
#' fitS <- coxph(Surv(time, event) ~ 1,
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' ## product limit
#' GS <- survfit(Surv(time, event) ~ 1, data = d)
#' range(predictCoxPL(fitS)$survival - GS$surv)
#'
#' fitS <- coxph(Surv(time, event) ~ strata(X2),
#'   data = d, ties = "breslow", x = TRUE, y = TRUE
#' )
#'
#' ## product limit
#' GS <- survfit(Surv(time, event) ~ X2, data = d)
#' range(predictCoxPL(fitS)$survival - GS$surv)
#'
## * predictCoxPL (code)
#' @rdname predictCoxPL
#' @export
predictCoxPL <- function(object,
                         times,
                         newdata = NULL,
                         type = c("cumhazard", "survival"),
                         keep.strata = TRUE,
                         keep.infoVar = FALSE,
                         ...) {
  stop <- NULL ## [:CRANcheck:] data.table

  ## ** normalize arguments
  object.modelFrame <- coxModelFrame(object)
  if (!is.null(newdata)) {
    if (is.data.table(newdata)) {
      newdata <- copy(newdata)
    } else {
      setDT(copy(newdata))
    }
  } else if (missing(times)) {
    times <- numeric(0)
  }

  if ("survival" %in% type == FALSE) {
    stop(
      "The argument \'type\' must contain \"survival\" \n",
      "use the function predictCox otherwise \n"
    )
  }

  ## ** run original prediction
  original.res <- predictCox(
    object = object,
    newdata = newdata,
    times = times,
    type = type,
    keep.strata = TRUE,
    keep.infoVar = TRUE,
    ...
  )

  infoVar <- original.res$infoVar
  if (keep.infoVar == FALSE) {
    original.res$infoVar <- NULL
  }

  ## no event: return exponential approximation (which equals 1 everywhere)
  if (all(object.modelFrame$status == 0)) {
    return(original.res)
  }

  ## ** compute survival
  if (is.null(newdata)) {
    hazard.res <- as.data.table(predictCox(
      object = object,
      type = "hazard",
      keep.strata = TRUE,
      ...
    ))

    if (!is.null(hazard.res$strata)) {
      Ustrata <- levels(hazard.res$strata)
      n.Ustrata <- length(Ustrata)

      for (iS in 1:n.Ustrata) { ## iS <- 1
        iIndex.all <- which(hazard.res$strata == Ustrata[iS])
        iIndex.times <- which(original.res$strata == Ustrata[iS])

        vec.surv <- c(1, cumprod(1 - hazard.res$hazard[iIndex.all]))

        index.jump <- prodlim::sindex(
          eval.times = original.res$times[iIndex.times],
          jump.times = c(0, hazard.res$times[iIndex.all])
        )

        original.res$survival[iIndex.times] <- vec.surv[index.jump]
      }
      if (keep.strata == FALSE) {
        original.res$strata <- NULL
      }
    } else {
      vec.surv <- c(1, cumprod(1 - hazard.res$hazard))

      index.jump <- prodlim::sindex(
        eval.times = original.res$times,
        jump.times = c(0, hazard.res$times)
      )

      original.res$survival <- vec.surv[index.jump]
    }

    if ("hazard" %in% type == FALSE) {
      original.res$hazard <- NULL
    }
  } else if (infoVar$is.strata) {
    object.strata <- coxStrata(object,
      data = NULL,
      strata.vars = infoVar$stratavars
    )
    new.strata <- original.res$strata
    if (keep.strata == FALSE) {
      original.res$strata <- NULL
    }
    Ustrata <- unique(new.strata)
    n.Ustrata <- length(Ustrata)

    for (iStrata in 1:n.Ustrata) { # iStrata <- 1
      indexStrata.object <- which(object.strata == Ustrata[iStrata])
      indexStrata.newdata <- which(new.strata == Ustrata[iStrata])

      all.times <- object.modelFrame[indexStrata.object, .SD$stop]
      all.times <- sort(unique(all.times[all.times <= max(times)]))

      if (length(all.times) > 0) {
        res.tempo <- predictCox(object,
          newdata = newdata[indexStrata.newdata, ],
          times = all.times,
          type = "hazard"
        )

        lastEventTime.tempo <- original.res$lastEventTime[levels(Ustrata) == Ustrata[iStrata]]
        if (0 %in% all.times) {
          hazard.tempo <- cbind(res.tempo$hazard, NA)
          all.timesA <- c(all.times, lastEventTime.tempo + 1e-10)
        } else {
          hazard.tempo <- cbind(0, res.tempo$hazard, NA)
          all.timesA <- c(0, all.times, lastEventTime.tempo + 1e-10)
        }
        if (original.res$diag) {
          index.jump <- prodlim::sindex(eval.times = times[indexStrata.newdata], jump.times = all.timesA)
          original.res$survival[indexStrata.newdata, ] <- cbind(sapply(seq_len(length(indexStrata.newdata)), function(iP) { # iP <- 16
            prod(1 - hazard.tempo[iP, 1:index.jump[iP]])
          }))
        } else {
          index.jump <- prodlim::sindex(eval.times = times, jump.times = all.timesA)
          survival.PL <- t(apply(hazard.tempo, 1, function(x) {
            cumprod(1 - x)
          }))
          original.res$survival[indexStrata.newdata, ] <- survival.PL[, index.jump, drop = FALSE]
        }
      }
    }
  } else {
    all.times <- unique(sort(object.modelFrame[stop <= max(times), stop]))
    if (length(all.times) > 0) {
      res.tempo <- predictCox(object,
        newdata = newdata,
        times = all.times,
        type = "hazard"
      )

      lastEventTime.tempo <- original.res$lastEventTime

      if (0 %in% all.times) {
        hazard.tempo <- cbind(res.tempo$hazard, NA)
        all.timesA <- c(all.times, lastEventTime.tempo + 1e-10)
      } else {
        hazard.tempo <- cbind(0, res.tempo$hazard, NA)
        all.timesA <- c(0, all.times, lastEventTime.tempo + 1e-10)
      }

      index.jump <- prodlim::sindex(eval.times = times, jump.times = all.timesA)
      if (original.res$diag) {
        original.res$survival <- cbind(sapply(seq_len(NROW(newdata)), function(iP) { # iP <- 16
          prod(1 - hazard.tempo[iP, 1:index.jump[iP]])
        }))
      } else {
        survival.PL <- t(apply(hazard.tempo, 1, function(x) {
          cumprod(1 - x)
        }))
        original.res$survival <- survival.PL[, index.jump, drop = FALSE]
      }
    }
  }

  ## ** export
  if (any(na.omit(as.double(original.res$survival)) > 1) || any(na.omit(as.double(original.res$survival)) < 0)) {
    warning(
      "Estimated survival outside the range [0,1]. \n",
      "Consider using predictCox instead of predictCoxPL. \n"
    )
  }
  return(original.res)
}
#----------------------------------------------------------------------
### predictCoxPL.R ends here
