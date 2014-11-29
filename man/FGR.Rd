\name{FGR}
\alias{FGR}
\title{Formula wrapper for crr from cmprsk}
\usage{
FGR(formula, data, cause = 1, y = TRUE, ...)
}
\arguments{
  \item{formula}{A formula whose left hand side is a
  \code{Hist} object -- see \code{\link{Hist}}.  The right
  hand side specifies (a linear combination of) the
  covariates. See examples below.}

  \item{data}{A data.frame in which all the variables of
  \code{formula} can be interpreted.}

  \item{cause}{The failure type of interest. Defaults to
  \code{1}.}

  \item{y}{logical value: if \code{TRUE}, the response
  vector is returned in component \code{response}.}

  \item{\dots}{...}
}
\value{
See \code{crr}.
}
\description{
Formula interface for Fine-Gray regression competing risk
models.
}
\details{
Formula interface for the function \code{crr} from the
\code{cmprsk} package.

The function \code{crr} allows to multiply some covariates
by time before they enter the linear predictor. This can be
achieved with the formula interface, however, the code
becomes a little cumbersome. See the examples.
}
\examples{
library(riskRegression)
library(prodlim)
library(pec)
library(cmprsk)
library(lava)
d <- SimCompRisk(100)
f1 <- FGR(Hist(time,cause)~X1+X2,data=d)
print(f1)

## crr allows that some covariates are multiplied by
## a function of time (see argument tf of crr)
## by FGR uses the identity matrix
f2 <- FGR(Hist(time,cause)~cov2(X1)+X2,data=d)
print(f2)

## same thing, but more explicit:
f3 <- FGR(Hist(time,cause)~cov2(X1)+cov1(X2),data=d)
print(f3)

## both variables can enter cov2:
f4 <- FGR(Hist(time,cause)~cov2(X1)+cov2(X2),data=d)
print(f4)

## change the function of time
qFun <- function(x){x^2}
noFun <- function(x){x}
sqFun <- function(x){x^0.5}

## multiply X1 by time^2 and X2 by time:
f5 <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2),data=d)

## the same but more explicit
f5a <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=noFun),data=d)

## multiply X1 by time^2 and X2 by sqrt(time)
f5b <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=sqFun),data=d,cause=1)

## additional arguments for crr
f6<- FGR(Hist(time,cause)~X1+X2,data=d, cause=1,gtol=1e-5)
}
\author{
Thomas Alexander Gerds \email{tag@biostat.ku.dk}
}
\references{
Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute
risk regression for competing risks: interpretation, link
functions and prediction Research report 11/7. Department
of Biostatistics, University of Copenhagen
}
\seealso{
\code{\link{riskRegression}}
}
\keyword{survival}

