% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.ate.R
\name{summary.ate}
\alias{summary.ate}
\title{Summary Average Treatment Effects}
\usage{
\method{summary}{ate}(
  object,
  estimator = object$estimator[1],
  short = FALSE,
  type = c("meanRisk", "diffRisk"),
  se = FALSE,
  quantile = FALSE,
  estimate.boot = TRUE,
  digits = 3,
  ...
)
}
\arguments{
\item{object}{object obtained with function \code{ate}}

\item{estimator}{[character] The type of estimator relative to which the estimates should be displayed.}

\item{short}{[logical] If \code{TRUE}, only displays the estimated risks.}

\item{type}{[character vector] what to displayed.
Can be \code{"meanRisk"} to display the risks specific to each treatment group,
\code{"diffRisk"} to display the difference in risks between treatment groups,
or \code{"ratioRisk"} to display the ratio of risks between treatment groups,.}

\item{se}{[logical] should the standard error of the risks be displayed?}

\item{quantile}{[logical] should the quantile of the confidence bands be displayed?}

\item{estimate.boot}{[logical] should the average estimate on the bootstrap samples be displayed?}

\item{digits}{[integer, >0] Number of digits.}

\item{...}{passed to confint}
}
\description{
Summary average treatment effects.
}
\details{
to display confidence intervals/bands and p.value,
the \code{confint} method needs to be applied on the object.
}
\seealso{
\code{\link[data.table]{as.data.table}} to extract the estimates in a \code{data.table} object.
\code{\link{autoplot.ate}} for a graphical representation the standardized risks.
\code{\link{confint.ate}} to compute p-values and adjusted p-values
or perform statistical inference using a transformation.
\code{\link{confint.ate}} to compute (pointwise/simultaneous) confidence intervals and (unadjusted/adjusted) p-values, possibly using a transformation.
}
