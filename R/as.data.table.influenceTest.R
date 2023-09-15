### as.data.table.influenceTest.R ---
## ----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  1 2018 (13:41)
## Version:
## Last-Updated: Jan 29 2019 (10:49)
##           By: Thomas Alexander Gerds
##     Update #: 22
## ----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
## ----------------------------------------------------------------------
##
### Code:

## * as.data.table.influenceTest (documentation)
#' @title Turn influenceTest Object Into a \code{data.table}
#' @description Turn influenceTest object into a \code{data.table}.
#' @name as.data.table.influenceTest
#'
#' @param x object obtained with function \code{influenceTest}
#' @param keep.rownames Not used.
#' @param se [logical] Should standard errors/quantile for confidence bands be displayed?
#' @param ... Not used.
#'

## * as.data.table.influenceTest (code)
#' @rdname as.data.table.influenceTest
#' @export
as.data.table.influenceTest <- function(x, keep.rownames = FALSE, se = TRUE, ...) {
  n.obs <- NROW(x$delta)
  n.time <- length(x$time)


  ls.newdata <- lapply(1:n.obs, function(iRow) {
    if (is.null(x$newdata)) {
      return(rep(iRow), n.time)
    } else {
      do.call(rbind, rep(list(x$newdata[iRow]), n.time))
    }
  })
  out <- data.table::rbindlist(ls.newdata)
  if (!is.null(x$strata)) {
    out$strata <- unlist(lapply(x$strata, rep, n.time))
  }
  out$time <- rep(x$time, n.obs)
  out$difference <- as.numeric(t(x$delta))

  if (se) {
    out$se <- as.numeric(t(x$delta.se))
  }
  if (!is.null(x$conf.level)) {
    out$lower <- as.numeric(t(x$delta.lower))
    out$upper <- as.numeric(t(x$delta.upper))
    out$p.value <- as.numeric(t(x$delta.p.value))
  }
  if (x$band[[1]] && !is.null(x$conf.level[[1]])) {
    if (se) {
      out$quantileBand <- unlist(lapply(x$delta.quantileBand, rep, n.time))
    }
    out$lowerBand <- as.numeric(t(x$delta.lowerBand))
    out$upperBand <- as.numeric(t(x$delta.upperBand))
  }
  return(out)
}


######################################################################
### as.data.table.influenceTest.R ends here
