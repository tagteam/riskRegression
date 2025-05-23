% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxBaseEstimator.R
\name{coxBaseEstimator}
\alias{coxBaseEstimator}
\alias{coxBaseEstimator.coxph}
\alias{coxBaseEstimator.phreg}
\alias{coxBaseEstimator.prodlim}
\alias{coxBaseEstimator.GLMnet}
\title{Extract the type of estimator for the baseline hazard}
\usage{
coxBaseEstimator(object)

\method{coxBaseEstimator}{coxph}(object)

\method{coxBaseEstimator}{phreg}(object)

\method{coxBaseEstimator}{prodlim}(object)

\method{coxBaseEstimator}{GLMnet}(object)
}
\arguments{
\item{object}{The fitted Cox regression model object either
obtained with \code{coxph} (survival package), \code{cph}
(rms package), or \code{phreg} (mets package).}
}
\description{
Extract the type of estimator for the baseline hazard
}
\author{
Brice Ozenne broz@sund.ku.dk
}
