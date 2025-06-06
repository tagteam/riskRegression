% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxVarCov.R
\name{coxVarCov}
\alias{coxVarCov}
\alias{coxVarCov.cph}
\alias{coxVarCov.coxph}
\alias{coxVarCov.phreg}
\alias{coxVarCov.prodlim}
\title{Extract the variance covariance matrix of the beta from a Cox model}
\usage{
coxVarCov(object)

\method{coxVarCov}{cph}(object)

\method{coxVarCov}{coxph}(object)

\method{coxVarCov}{phreg}(object)

\method{coxVarCov}{prodlim}(object)
}
\arguments{
\item{object}{The fitted Cox regression model object either
obtained with \code{coxph} (survival package), \code{cph}
(rms package), or \code{phreg} (mets package).}
}
\description{
Extract the variance covariance matrix of the beta from a Cox model
}
\details{
Should return \code{NULL} if the Cox model has no covariate. 
The rows and columns of the variance covariance matrix must be named with the names used in the design matrix.
}
\author{
Brice Ozenne broz@sund.ku.dk
}
