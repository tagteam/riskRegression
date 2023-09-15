### sim.synth.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  7 2022 (13:33)
## Version:
## Last-Updated: Jul  7 2022 (14:26)
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

#' @title Simulating from a synthesized object
#'
#' @description Simulating from a synthesized object
#' @param object generated with \code{synthesize}
#' @param n sample size
#' @param drop.latent if \code{TRUE} remove the latent event times from the resulting data set.
#' @param ... additional arguments passed on to \code{lava::sim}
#' @export simsynth
#' @examples
#' library(survival)
#' m <- synthesize(Surv(time, status) ~ sex + age + bili, data = pbc)
#' simsynth(m, 10, drop.latent = TRUE)
#'
simsynth <- function(object, n = 200, drop.latent = FALSE, ...) {
  lava.object <- object$lava.object
  res <- lava::sim(lava.object, n, ...)
  labels <- object$labels
  for (var in names(labels)) {
    res[[var]] <- factor(res[[var]])
    levels(res[[var]]) <- labels[[var]]
  }
  # remove variables that would not be in the original data set
  if (drop.latent) {
    # remove latent times
    if (length(lava.object$attributes$eventHistory$time$latentTimes) > 0) {
      res <- res[, -match(lava.object$attributes$eventHistory$time$latentTimes, names(res), nomatch = 0)]
    }
    # remove dummy variables
    categories <- object$categories
    for (c in categories) {
      res <- res[, -grep(c, names(res))[-1]]
    }
  }
  res
}


######################################################################
### sim.synth.R ends here
