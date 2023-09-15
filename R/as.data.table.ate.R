### as.data.table.ate.R ---
## ----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28)
## Version:
## Last-Updated: mar  8 2023 (09:58)
##           By: Brice Ozenne
##     Update #: 188
## ----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
## ----------------------------------------------------------------------
##
### Code:

## * as.data.table.ate (documentation)
#' @title Turn ate Object Into a \code{data.table}
#' @description Turn ate object into a \code{data.table}.
#' @name as.data.table.ate
#'
#' @param x object obtained with function \code{ate}
#' @param keep.rownames Not used.
#' @param estimator [character] The type of estimator relative to which the estimates should be output.
#' @param type [character vector] The type of risk to export.
#' Can be \code{"meanRisk"} to export the risks specific to each treatment group,
#' \code{"diffRisk"} to export the difference in risks between treatment groups,
#' or \code{"ratioRisk"} to export the ratio of risks between treatment groups.
#' @param ... Not used.
#'

## * as.data.table.ate (code)
#' @rdname as.data.table.ate
#' @export
as.data.table.ate <- function(x, keep.rownames = FALSE, estimator = x$estimator, type = c("meanRisk", "diffRisk", "ratioRisk"), ...) {
  estimator <- match.arg(estimator, choices = x$estimator, several.ok = TRUE)
  type <- match.arg(type, choices = c("meanRisk", "diffRisk", "ratioRisk"), several.ok = TRUE)

  if (!is.null(x$allContrasts)) {
    # allContrasts <- x$allContrasts
    contrasts <- x$contrasts
  } else {
    contrasts <- x$contrasts
    # allContrasts <- utils::combn(contrasts, m = 2)
  }

  ## ** meanRisk
  if ("meanRisk" %in% type) {
    ## iIndexRow <- which((x$meanRisk$estimator %in% estimator) * (x$meanRisk$treatment %in% contrasts) == 1)

    ## meanRisk <- x$meanRisk[iIndexRow]
    out1 <- cbind(
      type = "meanRisk",
      estimator = x$meanRisk$estimator,
      time = x$meanRisk$time,
      level = x$meanRisk$treatment,
      x$meanRisk[, .SD, .SDcols = setdiff(names(x$meanRisk), c("estimator", "time", "treatment"))]
    )
    if (x$inference$p.value && any(type %in% c("diffRisk", "ratioRisk"))) {
      out1$p.value <- as.numeric(NA)
      if (x$inference$band) {
        out1$adj.p.value <- as.numeric(NA)
      }
    }
  } else {
    out1 <- NULL
  }


  ## ** diffRisk
  if ("diffRisk" %in% type) {
    ## iIndexRow <- which((x$diffRisk$estimator %in% estimator) * (interaction(x$diffRisk$A, x$diffRisk$B) %in% interaction(allContrasts[1, ], allContrasts[2, ])) == 1)
    ## diffRisk <- x$diffRisk[iIndexRow]
    out2 <- cbind(
      type = "diffRisk",
      estimator = x$diffRisk$estimator,
      time = x$diffRisk$time,
      level = paste0(x$diffRisk$A, ".", x$diffRisk$B),
      x$diffRisk[, .SD, .SDcols = setdiff(names(x$diffRisk), c("estimator", "time", "A", "B", "estimate.A", "estimate.B"))]
    )
  } else {
    out2 <- NULL
  }
  ## ** ratioRisk
  if ("ratioRisk" %in% type) {
    ## iIndexRow <- which((x$ratioRisk$estimator %in% estimator) * (interaction(x$ratioRisk$A, x$ratioRisk$B) %in% interaction(allContrasts[1, ], allContrasts[2, ])) == 1)
    ## ratioRisk <- x$ratioRisk[iIndexRow]
    out3 <- cbind(
      type = "ratioRisk",
      estimator = x$ratioRisk$estimator,
      time = x$ratioRisk$time,
      level = paste0(x$ratioRisk$A, ".", x$ratioRisk$B),
      x$ratioRisk[, .SD, .SDcols = setdiff(names(x$ratioRisk), c("estimator", "time", "A", "B", "estimate.A", "estimate.B"))]
    )
  } else {
    out3 <- NULL
  }

  ## ** export
  out <- rbind(out1, out2, out3)
  if (all(is.na(out$time)) && is.na(x$variable["time"])) {
    out[, c("time") := NULL]
  }
  return(out[])
}



######################################################################
### as.data.table.ate.R ends here
