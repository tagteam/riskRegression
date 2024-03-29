% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoplot.predictCox.R
\name{autoplot.predictCox}
\alias{autoplot.predictCox}
\title{Plot Predictions From a Cox Model}
\usage{
\method{autoplot}{predictCox}(
  object,
  type = NULL,
  ci = object$se,
  band = object$band,
  plot = TRUE,
  smooth = NULL,
  digits = 2,
  alpha = NA,
  group.by = "row",
  reduce.data = FALSE,
  ylab = NULL,
  first.derivative = FALSE,
  ...
)
}
\arguments{
\item{object}{Object obtained with the function \code{predictCox}.}

\item{type}{[character] The type of predicted value to display.
Choices are:
\code{"hazard"} the hazard function,
\code{"cumhazard"} the cumulative hazard function, 
or \code{"survival"} the survival function.}

\item{ci}{[logical] If \code{TRUE} display the confidence intervals for the predictions.}

\item{band}{[logical] If \code{TRUE} display the confidence bands for the predictions.}

\item{plot}{[logical] Should the graphic be plotted.}

\item{smooth}{[logical] Should a smooth version of the risk function be plotted instead of a simple function?}

\item{digits}{[integer] Number of decimal places when displaying the values of the covariates in the caption.}

\item{alpha}{[numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.}

\item{group.by}{[character] The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.}

\item{reduce.data}{[logical] If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.}

\item{ylab}{[character] Label for the y axis.}

\item{first.derivative}{[logical] If \code{TRUE}, display the first derivative over time of the risks/risk differences/risk ratios.
(confidence intervals are obtained via simulation).}

\item{...}{Additional parameters to cutomize the display.}
}
\value{
Invisible. A list containing:
\itemize{
\item plot: the ggplot object.
\item data: the data used to create the plot.
}
}
\description{
Plot predictions from a Cox model.
}
\examples{
library(survival)
library(ggplot2)

#### simulate data ####
set.seed(10)
d <- sampleData(1e2, outcome = "survival")
seqTau <- c(0,sort(unique(d$time[d$event==1])), max(d$time))

#### Cox model ####
m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
                data = d, x = TRUE, y = TRUE)

## display baseline hazard
e.basehaz <- predictCox(m.cox)
autoplot(e.basehaz, type = "cumhazard")
\dontrun{
autoplot(e.basehaz, type = "cumhazard", size.point = 0) ## without points
autoplot(e.basehaz, type = "cumhazard", smooth = TRUE)
autoplot(e.basehaz, type = "cumhazard", smooth = TRUE, first.derivative = TRUE)
}

## display baseline hazard with type of event
\dontrun{
e.basehaz <- predictCox(m.cox, keep.newdata = TRUE)
autoplot(e.basehaz, type = "cumhazard")
autoplot(e.basehaz, type = "cumhazard", shape.point = c(3,NA))
}

## display predicted survival
\dontrun{
pred.cox <- predictCox(m.cox, newdata = d[1:2,],
  times = seqTau, type = "survival", keep.newdata = TRUE)
autoplot(pred.cox)
autoplot(pred.cox, smooth = TRUE)
autoplot(pred.cox, group.by = "covariates")
autoplot(pred.cox, group.by = "covariates", reduce.data = TRUE)
autoplot(pred.cox, group.by = "X1", reduce.data = TRUE)
}

## predictions with confidence interval/bands
\dontrun{
pred.cox <- predictCox(m.cox, newdata = d[1:2,,drop=FALSE],
  times = seqTau, type = "survival", band = TRUE, se = TRUE, keep.newdata = TRUE)
res <- autoplot(pred.cox, ci = TRUE, band = TRUE, plot = FALSE)
res$plot + facet_wrap(~row)
res2 <- autoplot(pred.cox, ci = TRUE, band = TRUE, alpha = 0.1, plot = FALSE)
res2$plot + facet_wrap(~row)
}

#### Stratified Cox model ####
\dontrun{
m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X4,
                      data = d, x = TRUE, y = TRUE)

## baseline hazard
pred.baseline <- predictCox(m.cox.strata, keep.newdata = TRUE, type = "survival")
res <- autoplot(pred.baseline)
res$plot + facet_wrap(~strata, labeller = label_both)

## predictions
pred.cox.strata <- predictCox(m.cox.strata, newdata = d[1:3,,drop=FALSE],
                              time = seqTau, keep.newdata = TRUE, se = TRUE)

res2 <- autoplot(pred.cox.strata, type = "survival", group.by = "strata", plot = FALSE)
res2$plot + facet_wrap(~strata, labeller = label_both) + theme(legend.position="bottom")

## smooth version
autoplot(pred.cox.strata, type = "survival", group.by = "strata", smooth = TRUE, ci = FALSE)
}

#### Cox model with splines ####
\dontrun{
require(splines)
m.cox.spline <- coxph(Surv(time,event)~ X1 + X2 + ns(X6,4),
                data = d, x = TRUE, y = TRUE)
grid <- data.frame(X1 = factor(0,0:1), X2 = factor(0,0:1),
                   X6 = seq(min(d$X6),max(d$X6), length.out = 100))
pred.spline <- predictCox(m.cox.spline, newdata = grid, keep.newdata = TRUE,
                          se = TRUE, band = TRUE, centered = TRUE, type = "lp")
autoplot(pred.spline, group.by = "X6")
autoplot(pred.spline, group.by = "X6", alpha = 0.5)

grid2 <- data.frame(X1 = factor(1,0:1), X2 = factor(0,0:1),
                    X6 = seq(min(d$X6),max(d$X6), length.out = 100))
pred.spline <- predictCox(m.cox.spline, newdata = rbind(grid,grid2), keep.newdata = TRUE,
                          se = TRUE, band = TRUE, centered = TRUE, type = "lp")
autoplot(pred.spline, group.by = c("X6","X1"), alpha = 0.5, plot = FALSE)$plot + facet_wrap(~X1)
}
}
\seealso{
\code{\link{predictCox}} to compute cumulative hazard and survival based on a Cox model.
}
