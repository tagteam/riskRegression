\name{subjectWeights}
\alias{subjectWeights}
\alias{subjectWeights.aalen}
\alias{subjectWeights.cox}
\alias{subjectWeights.km}
\alias{subjectWeights.marginal}
\alias{subjectWeights.none}
\alias{subjectWeights.nonpar}
\title{Estimation of censoring probabilities at subject specific times}
\usage{
subjectWeights(formula, data, method = c("cox", "marginal", "km", "nonpar",
  "forest", "none"), args, lag = 1)
}
\arguments{
  \item{formula}{A survival formula like,
  Surv(time,status)~1 or Hist(time,status)~1 where status=0
  means censored. The status variable is internally
  reversed for estimation of censoring rather than survival
  probabilities. Some of the available models, see argument
  \code{model}, will use predictors on the right hand side
  of the formula.}

  \item{data}{The data used for fitting the censoring
  model}

  \item{method}{Censoring model used for estimation of the
  (conditional) censoring distribution.}

  \item{args}{Arguments passed to the fitter of the
  method.}

  \item{lag}{If equal to \code{1} then obtain
  \code{G(T_i-|X_i)}, if equal to \code{0} estimate the
  conditional censoring distribution at the subjectTimes,
  i.e. (\code{G(T_i|X_i)}).}
}
\value{
\item{times}{The times at which weights are estimated}
\item{weights}{Estimated weights at individual time values
\code{subjectTimes}} \item{lag}{The time lag.}
\item{fit}{The fitted censoring model} \item{method}{The
method for modelling the censoring distribution}
\item{call}{The call}
}
\description{
This function is used internally to contruct pseudo values
by inverse of the probability of censoring weights.
}
\details{
Inverse of the probability of censoring weights usually
refer to the probabilities of not being censored at certain
time points. These probabilities are also the values of the
conditional survival function of the censoring time given
covariates. The function subjectWeights estimates the
conditional survival function of the censoring times and
derives the weights.

IMPORTANT: the data set should be ordered,
\code{order(time,-status)} in order to get the
\code{weights} in the right order for some choices of
\code{method}.
}
\examples{
library(pec)
library(prodlim)
library(survival)
dat=SimSurv(300)

dat <- dat[order(dat$time,-dat$status),]

# using the marginal Kaplan-Meier for the censoring times

WKM=subjectWeights(Hist(time,status)~X2,data=dat,method="marginal")
plot(WKM$fit)
WKM$fit
WKM$weights

# using the Cox model for the censoring times given X2

WCox=subjectWeights(Surv(time,status)~X2,data=dat,method="cox")
WCox
plot(WCox$weights,WKM$weights)

# using the stratified Kaplan-Meier for the censoring times given X2

WKM2 <- subjectWeights(Surv(time,status)~X2,data=dat,method="nonpar")
plot(WKM2$fit,add=FALSE)
}
\author{
Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\keyword{survival}

