### plotAUC.R ---
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (09:19)
## Version:
## last-updated: May 30 2023 (08:07)
##           By: Thomas Alexander Gerds
##     Update #: 86
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Plot of time-dependent AUC curves
##'
##' @title Plot of time-dependent AUC curves
##' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
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
##' @examples
##' set.seed(9)
##' library(survival)
##' library(prodlim)
##' set.seed(10)
##' d=sampleData(100,outcome="survival")
##' nd=sampleData(100,outcome="survival")
##' f1=coxph(Surv(time,event)~X1+X6+X8,data=d,x=TRUE,y=TRUE)
##' f2=coxph(Surv(time,event)~X2+X5+X9,data=d,x=TRUE,y=TRUE)
##' xx=Score(list("X1+X6+X8"=f1,"X2+X5+X9"=f2), formula=Surv(time,event)~1,
##' data=nd, metrics="auc", null.model=FALSE, times=seq(3:10))
##' aucgraph <- plotAUC(xx)
##' plotAUC(xx,conf.int=TRUE)
##' ## difference between
##' plotAUC(xx,which="contrasts",conf.int=TRUE)
##'
#'
#' @export
plotAUC <- function(x,
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
  times <- contrast <- model <- AUC <- lower <- upper <- lower <- upper <- delta.AUC <- reference <- se <- NULL
  ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  pframe <- switch(which,
    "score" = {
      copy(x$AUC$score)
    },
    "contrasts" = {
      copy(x$AUC$contrasts)
    },
    {
      stop("argument 'which' has to be either 'score' for AUC or 'contrasts' for differences in AUC.")
    }
  )
  if (length(pframe$times) < 2) stop(paste("Need at least two time points for plotting time-dependent AUC. Object has only ", length(pframe$times), "times"))
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
    pframe[is.na(se) & times == 0, lower := 0]
    pframe[is.na(se) & times == 0, upper := 0]
  } else {
    pframe[, contrast := factor(paste(model, reference, sep = " - "))]
    mm <- unique(pframe$contrast)
    pframe[is.na(se) & times == 0, lower := 0]
    pframe[is.na(se) & times == 0, upper := 0]
  }
  lenmm <- length(mm)
  if (missing(xlab)) xlab <- "Time"
  if (missing(ylab)) if (which == "score") ylab <- "AUC" else ylab <- expression(paste(Delta, " AUC"))
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
    if (which == "score") {
      ylim <- c(0.5, 1)
      axis2.DefaultArgs <- list(side = 2, las = 2, at = seq(0, ylim[2], ylim[2] / 4), mgp = c(4, 1, 0))
    } else {
      if (is.null(pframe$lower) || all(is.na(pframe$lower))) {
        ylim <- c(-1, 1)
      } else {
        ylim <- c(floor(10 * min(pframe$lower)) / 10, ceiling(10 * max(pframe$upper)) / 10)
      }
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
  if (which == "score") {
    legend.DefaultArgs <- list(legend = mm, lwd = lwd, col = col, lty = lty, cex = cex, bty = "n", y.intersp = 1.3, x = "topleft")
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
  if (which == "score") {
    ## AUC
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
      thisline$y <- pframe[model == mod][["AUC"]]
      do.call("lines", thisline)
    }
  } else {
    ## delta AUC
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
      thisline$y <- pframe[contrast == con][["delta.AUC"]]
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
### plotAUC.R ends here
