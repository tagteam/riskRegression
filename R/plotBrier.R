### plotBrier.R ---
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 23 2017 (11:07)
## Version:
## last-updated: May 30 2023 (08:07)
##           By: Thomas Alexander Gerds
##     Update #: 76
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Plot Brier score curves
##' @title Plot Brier curve
#' @param x Object obtained with \code{Score}
#' @param models Choice of models to plot
#' @param which Character. Either \code{"score"} to show AUC or
#'     \code{"contrasts"} to show differences between AUC.
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param col line color
#' @param lwd line width
#' @param lty line style
#' @param cex point size
#' @param pch point style
#' @param type line type
#' @param axes Logical. If \code{TRUE} draw axes.
#' @param percent Logical. If \code{TRUE} scale y-axis in percent.
#' @param conf.int Logical. If \code{TRUE} draw confidence shadows.
#' @param legend Logical. If \code{TRUE} draw legend.
#' @param ... Used for additional control of the subroutines: plot,
#'     axis, lines, legend. See \code{\link{SmartControl}}.
##' @examples
##' # survival
##' library(survival)
##' library(prodlim)
##' set.seed(7)
##' ds1=sampleData(40,outcome="survival")
##' ds2=sampleData(40,outcome="survival")
##' f1 <- coxph(Surv(time,event)~X1+X3+X5+X7+X9,data=ds1,x=TRUE)
##' f2 <- coxph(Surv(time,event)~X2+X4+X6+X8+X10,data=ds1,x=TRUE)
##' xscore <- Score(list(f1,f2),formula=Hist(time,event)~1,data=ds2,times=0:12,metrics="brier")
##' plotBrier(xscore)
#' @export
#'
#'
plotBrier <- function(x,
                      models,
                      which = "score",
                      xlim,
                      ylim,
                      xlab,
                      ylab,
                      col,
                      lwd,
                      lty = 1,
                      cex = 1,
                      pch = 1,
                      type = "l",
                      axes = 1L,
                      percent = 1L,
                      conf.int = 0L,
                      legend = 1L,
                      ...) {
  times <- contrast <- model <- se <- lower <- upper <- reference <- NULL
  ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  which <- tolower(which[1])
  pframe <- switch(which,
    "score" = {
      copy(x$Brier$score)
    },
    "ipa" = {
      copy(x$Brier$score)
    },
    "contrasts" = {
      copy(x$Brier$contrasts)
    },
    {
      stop("argument 'which' has to be either 'score' for Brier, 'ipa' for IPA, or 'contrasts' for differences in Brier.")
    }
  )
  if (length(pframe$times) < 2) stop(paste("Need at least two time points for plotting time-dependent Brier. Object has only ", length(pframe$times), "times"))
  if (!missing(models)) {
    if (any(!(models %in% unique(pframe$models)))) {
      stop(paste0(
        "Fitted object does not contain models named: ", paste0(models, collapse = ", "),
        "\nAvailable models are named: ", paste0(unique(pframe$models), collapse = ", ")
      ))
    }
    pframe <- pframe[model %in% models]
  }
  if (which == "score") {
    mm <- unique(pframe$model)
    if ("se" %in% names(pframe)) {
      pframe[is.na(se) & times == 0, lower := 0]
      pframe[is.na(se) & times == 0, upper := 0]
    }
  } else {
    if (which == "ipa") {
      mm <- unique(pframe[model != "Null model", model[1]])
      pframe <- pframe[model %in% mm]
      if ("se" %in% names(pframe)) {
        pframe[is.na(se) & times == 0, lower := 0]
        pframe[is.na(se) & times == 0, upper := 0]
      }
    } else {
      pframe[, contrast := factor(paste(model, reference, sep = " - "))]
      mm <- unique(pframe$contrast)
      if ("se" %in% names(pframe)) {
        pframe[is.na(se) & times == 0, lower := 0]
        pframe[is.na(se) & times == 0, upper := 0]
      }
    }
  }
  lenmm <- length(mm)
  if (missing(xlab)) xlab <- "Time"
  if (missing(ylab)) {
    ylab <- switch(which,
      "score" = "Brier score",
      "ipa" = "Index of prediction accuracy",
      expression(paste(Delta, " Brier score"))
    )
  }
  if (missing(col)) col <- rep(cbbPalette, length.out = lenmm)
  names(col) <- mm
  if (missing(lwd)) lwd <- 2
  lwd <- rep(lwd, length.out = lenmm)
  names(lwd) <- mm
  pch <- rep(pch, length.out = lenmm)
  names(pch) <- mm
  type <- rep(type, length.out = lenmm)
  names(type) <- mm
  if (missing(lwd)) lty <- 1
  lty <- rep(lty, length.out = lenmm)
  names(lty) <- mm
  if (missing(xlim)) xlim <- pframe[, range(times)]
  if (missing(ylim)) {
    if (which %in% c("score", "ipa")) {
      ylim <- c(0, .3)
      if (which == "ipa") {
        ylim <- c(0, max(pframe$IPA, na.rm = 0))
      }
      axis2.DefaultArgs <- list(
        side = 2,
        las = 2,
        at = seq(0, ylim[2], ylim[2] / 4),
        mgp = c(4, 1, 0)
      )
    } else {
      ylim <- c(floor(10 * min(pframe$lower)) / 10, ceiling(10 * max(pframe$upper)) / 10)
      yat <- seq(ylim[1], ylim[2], 0.05)
      ## this is a strange behaviour of R: seq(-0.6,.1,0.05)
      ## [1] -6.000000e-01 -5.500000e-01 -5.000000e-01 -4.500000e-01 -4.000000e-01 -3.500000e-01 -3.000000e-01 -2.500000e-01
      ## [9] -2.000000e-01 -1.500000e-01 -1.000000e-01 -5.000000e-02  1.110223e-16  5.000000e-02  1.000000e-01
      yat <- round(100 * yat) / 100
      ## axis2.DefaultArgs <- list(side=2,las=2,at=seq(ylim[1],ylim[2],abs(ylim[2]-ylim[1])/4),mgp=c(4,1,0))
      axis2.DefaultArgs <- list(side = 2, las = 2, at = yat, mgp = c(4, 1, 0))
    }
  } else {
    axis2.DefaultArgs <- list(side = 2, las = 2, at = seq(ylim[1], ylim[2], abs(ylim[2] - ylim[1]) / 4), mgp = c(4, 1, 0))
  }
  lines.DefaultArgs <- list(pch = pch, type = type, cex = cex, lwd = lwd, col = col, lty = lty)
  axis1.DefaultArgs <- list(side = 1, las = 1, at = seq(0, xlim[2], xlim[2] / 4))
  if (which %in% c("score", "ipa")) {
    legend.DefaultArgs <- list(legend = mm, lwd = lwd, col = col, lty = lty, cex = cex, bty = "n", y.intersp = 1.3, x = "topleft")
    if (which == "ipa") legend.DefaultArgs$x <- "bottomleft"
  } else {
    legend.DefaultArgs <- list(legend = as.character(unique(pframe$contrast)), lwd = lwd, col = col, lty = lty, cex = cex, bty = "n", y.intersp = 1.3, x = "topleft")
  }
  plot.DefaultArgs <- list(x = 0, y = 0, type = "n", ylim = ylim, xlim = xlim, ylab = ylab, xlab = xlab)
  control <- prodlim::SmartControl(
    call = list(...),
    keys = c("plot", "lines", "legend", "axis1", "axis2"),
    ignore = NULL,
    ignore.case = TRUE,
    defaults = list(
      "plot" = plot.DefaultArgs,
      "lines" = lines.DefaultArgs,
      "legend" = legend.DefaultArgs,
      "axis1" = axis1.DefaultArgs,
      "axis2" = axis2.DefaultArgs
    ),
    forced = list(
      "plot" = list(axes = FALSE),
      "axis1" = list(side = 1)
    ),
    verbose = TRUE
  )

  do.call("plot", control$plot)
  if (which %in% c("score", "ipa")) {
    ## Brier
    # not a very nice solution but this fixes the problem with the plot
    model.list <- unique(pframe[["model"]])
    times <- unique(pframe[["times"]])
    for (mod in model.list) {
      thisline <- control$line
      thisline$col <- thisline$col[[as.character(mod)]]
      thisline$lwd <- thisline$lwd[[as.character(mod)]]
      thisline$lty <- thisline$lty[[as.character(mod)]]
      thisline$pch <- thisline$pch[[as.character(mod)]]
      thisline$type <- thisline$type[[as.character(mod)]]
      thisline$x <- times
      if (which == "ipa") {
        thisline$y <- pframe[model == mod][["IPA"]] ## does this work?
      } else {
        thisline$y <- pframe[model == mod][["Brier"]]
      }
      do.call("lines", thisline)
    }
  } else {
    ## delta Brier
    # not a very nice solution but this fixes the problem with the plot
    contrast.list <- unique(pframe[["contrast"]])
    times <- unique(pframe[["times"]])
    for (con in contrast.list) {
      thisline <- control$line
      thisline$col <- thisline$col[[as.character(con)]]
      thisline$lwd <- thisline$lwd[[as.character(con)]]
      thisline$lty <- thisline$lty[[as.character(con)]]
      thisline$pch <- thisline$pch[[as.character(con)]]
      thisline$type <- thisline$type[[as.character(con)]]
      thisline$x <- times
      thisline$y <- pframe[contrast == con][["delta.Brier"]]
      do.call("lines", thisline)
    }
  }
  ## legend
  if (!(is.logical(legend[[1]]) && legend[[1]] == FALSE)) {
    do.call("legend", control$legend)
  }
  ## x-axis
  if (conf.int == TRUE) {
    dimcol <- sapply(col, function(cc) {
      prodlim::dimColor(cc)
    })
    names(dimcol) <- names(col)
    if (which == "score") {
      pframe[, polygon(x = c(times, rev(times)), y = c(lower, rev(upper)), col = dimcol[[as.character(model)]], border = NA), by = model]
    } else {
      pframe[, polygon(x = c(times, rev(times)), y = c(lower, rev(upper)), col = dimcol[[as.character(contrast)]], border = NA), by = contrast]
    }
  }
  if (axes) {
    control$axis2$labels <- paste(100 * control$axis2$at, "%")
    do.call("axis", control$axis1)
    do.call("axis", control$axis2)
  }
  invisible(pframe)
}


#----------------------------------------------------------------------
### plotBrier.R ends here
