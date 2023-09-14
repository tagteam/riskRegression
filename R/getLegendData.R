### getLegendData.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 13 2017 (10:50)
## Version:
## Last-Updated: Sep 19 2019 (12:27)
##           By: Thomas Alexander Gerds
##     Update #: 62
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
getLegendData <- function(object,
                          models,
                          times,
                          brier.in.legend = TRUE,
                          format.brier,
                          auc.in.legend = TRUE,
                          format.auc,
                          drop.null.model = TRUE,
                          scale = 100,
                          digits = 1,
                          ...) {
  model <- AUC <- lower <- upper <- Brier <- NULL
  if (missing(models)) {
    models <- names(object$models)
  }
  if (!is.null(object$null.model)) {
    if (drop.null.model == TRUE) models <- models[models != object$null.model]
  }
  maxlen <- max(nchar(as.character(models)))
  legend.text.models <- sprintf(paste0("%", maxlen, "s"), models)
  if (missing(format.auc)) {
    format.auc <- paste0("%1.", digits, "f [%1.", digits, "f;%1.", digits, "f]")
  }
  if (missing(format.brier)) {
    format.brier <- paste0("%1.", digits, "f [%1.", digits, "f;%1.", digits, "f]")
  }
  if (is.null(object$null.model) || drop.null.model[[1]] == TRUE) {
    keep.null.model <- FALSE
  } else {
    if (brier.in.legend == TRUE) {
      keep.null.model <- TRUE
    } else {
      keep.null.model <- match(object$null.model, models, nomatch = FALSE)
    }
  }
  if (auc.in.legend == TRUE) {
    if (is.null(object$AUC)) {
      warning("Cannot show AUC as it is not stored in object. Set metrics='auc' in the call of Score.")
      legend.text.auc <- NULL
    } else {
      auc.data <- object$AUC$score[(model %in% models)]
      if (object$response.type != "binary") {
        if (missing(times)) {
          tp <- max(auc.data[["times"]])
          if (length(unique(auc.data$times)) > 1) {
            warning("Time point not specified, use max of the available times: ", tp)
          }
        } else { ## can only do one time point
          tp <- times[[1]]
          if (!(tp %in% unique(auc.data$times))) {
            stop(paste0("Requested time ", times[[1]], " is not in object"))
          }
        }
        auc.data <- auc.data[times == tp]
      } else {
        tp <- NULL
      }
      if (keep.null.model == FALSE) {
        if (!is.null(object$null.model)) {
          auc.data <- auc.data[model != object$null.model]
        }
      }
      ## user's order
      auc.data[, model := factor(model, levels = models)]
      setkey(auc.data, model)
      legend.text.auc <- auc.data[, sprintf(fmt = format.auc, scale * AUC, scale * lower, scale * upper)]
    }
  } else {
    legend.text.auc <- NULL
  }
  if (brier.in.legend == TRUE) {
    if (is.null(object$Brier)) {
      warning("Cannot show Brier score as it is not stored in object. Set metrics='brier' in the call of Score.")
      legend.text.brier <- NULL
    } else {
      brier.data <- object$Brier$score[(model %in% models)]
      if (object$response.type != "binary") {
        if (missing(times)) {
          tp <- max(brier.data[["times"]])
          if (length(unique(brier.data$times)) > 1) {
            warning("Time point not specified, use max of the available times: ", tp)
          }
        } else { ## can only do one time point
          tp <- times[[1]]
          if (!(tp %in% unique(brier.data$times))) {
            stop(paste0("Requested time ", times[[1]], " is not in object"))
          }
        }
        brier.data <- brier.data[times == tp]
      } else {
        tp <- NULL
      }
      if (!is.null(object$null.model) && keep.null.model[[1]] == FALSE) {
        brier.data <- brier.data[model != object$null.model]
      }
      brier.data[, model := factor(model, levels = models)]
      setkey(brier.data, model)
      legend.text.brier <- brier.data[, sprintf(fmt = format.brier, scale * Brier, scale * lower, scale * upper)]
    }
  } else {
    legend.text.brier <- NULL
  }
  out <- cbind(legend.text.models,
    AUC = legend.text.auc,
    Brier = legend.text.brier
  )
  out
}


######################################################################
### getLegendData.R ends here
