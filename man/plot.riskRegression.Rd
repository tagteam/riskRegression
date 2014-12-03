\name{plot.riskRegression}
\alias{plot.CauseSpecificCox}
\alias{plot.predictedRisk}
\alias{plot.riskRegression}
\title{Plotting predicted risk}
\usage{
\method{plot}{riskRegression}(x, cause, newdata, xlab, ylab, xlim, ylim, lwd,
  col, lty, axes = TRUE, percent = TRUE, legend = TRUE, add = FALSE,
  ...)
}
\arguments{
  \item{x}{Fitted object obtained with one of \code{ARR},
  \code{LRR}, \code{riskRegression}.}

  \item{cause}{For CauseSpecificCox models the cause of
  interest.}

  \item{newdata}{A data frame containing predictor variable
  combinations for which to compute predicted risk.}

  \item{xlim}{See \code{plot}}

  \item{ylim}{See \code{plot}}

  \item{xlab}{See \code{plot}}

  \item{ylab}{See \code{plot}}

  \item{lwd}{A vector of line thicknesses for the
  regression coefficients.}

  \item{col}{A vector of colors for the regression
  coefficients.}

  \item{lty}{A vector of line types for the regression
  coefficients.}

  \item{axes}{Logical. If \code{FALSE} then do not draw
  axes.}

  \item{percent}{If true the y-axis is labeled in percent.}

  \item{legend}{If true draw a legend.}

  \item{add}{Logical. If \code{TRUE} then add lines to an
  existing plot.}

  \item{\dots}{Used for transclusion of smart arguments for
  \code{plot}, \code{lines}, \code{axis} and
  \code{background}. See function
  \code{\link{SmartControl}} from prodlim.}
}
\description{
Show predicted risk obtained by a risk prediction model as
a function of time.
}
\examples{
library(pec)
library(prodlim)
data(Melanoma)

fit.arr <- ARR(Hist(time,status)~invasion+age+strata(sex),data=Melanoma,cause=1)
plot(fit.arr)

fit.csc <- CSC(Hist(time,status)~invasion+age+sex,data=Melanoma,cause=1)
plot(fit.csc)

#'
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\keyword{survival}

