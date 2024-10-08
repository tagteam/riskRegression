% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictCoxPL.R
\name{predictCoxPL}
\alias{predictCoxPL}
\title{Deprecated Function for Product Limit Estimation of Survival Probabilities .}
\usage{
predictCoxPL(object, ...)
}
\arguments{
\item{object}{The fitted Cox regression model object either
obtained with \code{coxph} (survival package) or \code{cph}
(rms package).}

\item{...}{additional arguments to be passed to \code{\link{predictCox}}.}
}
\description{
Depreciated function for Product Limit Estimation of Survival Probabilities from a Cox model.
Use the \code{\link{predictCox}} function instead with argument \code{product.limit=TRUE}.
}
