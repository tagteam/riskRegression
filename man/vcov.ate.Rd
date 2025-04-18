% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcov.ate.R
\name{vcov.ate}
\alias{vcov.ate}
\title{Variance-Covariance Matrix for the Average Treatment Effect.}
\usage{
\method{vcov}{ate}(
  object,
  contrasts = NULL,
  times = NULL,
  estimator = NULL,
  type = NULL,
  ...
)
}
\arguments{
\item{object}{A \code{ate} object, i.e. output of the \code{ate} function.}

\item{contrasts}{[character vector] levels of the treatment variable for which the variance-covariance matrix should be assessed. Default is to consider all levels.}

\item{times}{[numeric vector] The timepoints at which the variance-covariance matrix should be displayed. Default is to consider all timepoints.}

\item{estimator}{[character] The type of estimator relative to which the variance-covariance matrix should be displayed.}

\item{type}{[character] should the variance-covariance matrix w.r.t. the average risk per treatment be displayed (\code{"meanRisk"}),
or the difference in average risk between any two pairs of treatments (\code{"diffRisk"}),
or the ratio in average risk between any two pairs of treatments (\code{"ratioRisk"}).}

\item{...}{Not used. For compatibility with the generic method.}
}
\value{
A numeric matrix.
}
\description{
Variance covariance matrix for the estimated average treatment effect.
}
\author{
Brice Ozenne \email{broz@sund.ku.dk}
}
