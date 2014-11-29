\name{riskRegression}
\alias{ARR}
\alias{LRR}
\alias{riskRegression}
\title{Risk Regression}
\usage{
riskRegression(formula, data, times, link = "relative", cause,
  confint = TRUE, cens.model, cens.formula, numSimu = 0, maxiter = 50,
  silent = 1, convLevel = 6, conservative = TRUE, ...)
}
\arguments{
  \item{formula}{Formula where the left hand side specifies
  the event history response and the right hand side the
  linear predictor.  See examples.}

  \item{data}{The data for fitting the model in which
  includes all the variables included in formula.}

  \item{times}{Vector of times. For each time point in
  \code{times} estimate the baseline risk and the
  timevarying coefficients.}

  \item{link}{\code{"relative"} for the absolute risk
  regression model.  \code{"logistic"} for the logistic
  risk regression model. \code{"prop"} for the Fine-Gray
  regression model.}

  \item{cause}{The cause of interest.}

  \item{confint}{If \code{TRUE} return the iid
  decomposition, that can be used to construct confidence
  bands for predictions.}

  \item{cens.model}{Specified the model for the
  (conditional) censoring distribution used for deriving
  weights (ICPW). Defaults to "KM" (the Kaplan-Meier method
  ignoring covariates) alternatively it may be "Cox" (Cox
  regression).}

  \item{cens.formula}{Right hand side of the formula used
  for fitting the censoring model.  If not specified the
  right hand side of \code{formula} is used.}

  \item{numSimu}{Number of simulations in resampling.}

  \item{maxiter}{Maximal number of iterations.}

  \item{silent}{Set this to \code{0} to see some system
  messages during fitting.}

  \item{convLevel}{Integer between 1 and 10, specifying the
  convergence criterion for the fitting algorithm as
  \code{10^{-convLevel}}.}

  \item{conservative}{If \code{TRUE} use variance formula
  that ignores the contribution by the estimate of the
  inverse of the probability of censoring weights}

  \item{...}{}

  \item{\dots}{Not used.}
}
\description{
Fits a regression model for the risk of an event --
allowing for competing risks.
}
\details{
This is more a twin-sister than a wrapper of the function
\code{\link{comp.risk}} (timereg package).
}
\examples{
data(Melanoma,package="riskRegression")
## tumor thickness on the log-scale
Melanoma$logthick <- log(Melanoma$thick)

# {{{ Single binary factor

## absolute risk regression
library(prodlim)
fit.arr <- ARR(Hist(time,status)~sex,data=Melanoma,cause=1)
print(fit.arr)
# show predicted cumulative incidences
plot(fit.arr,col=3:4,newdata=data.frame(sex=c("Female","Male")))

## compare with non-parametric Aalen-Johansen estimate
library(prodlim)
fit.aj <- prodlim(Hist(time,status)~sex,data=Melanoma)
plot(fit.aj,confint=FALSE)
plot(fit.arr,add=TRUE,col=3:4,newdata=data.frame(sex=c("Female","Male")))

## with time-dependent effect
fit.tarr <- ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1)
plot(fit.tarr,newdata=data.frame(sex=c("Female","Male")))

## logistic risk regression
fit.lrr <- LRR(Hist(time,status)~sex,data=Melanoma,cause=1)
summary(fit.lrr)

# }}}

# {{{ Single continuous factor

## tumor thickness on the log-scale
Melanoma$logthick <- log(Melanoma$thick)

## absolute risk regression
fit2.arr <- ARR(Hist(time,status)~logthick,data=Melanoma,cause=1)
print(fit2.arr)
# show predicted cumulative incidences
plot(fit2.arr,col=1:5,newdata=data.frame(logthick=quantile(Melanoma$logthick)))

## comparison with nearest neighbor non-parametric Aalen-Johansen estimate
library(prodlim)
fit2.aj <- prodlim(Hist(time,status)~logthick,data=Melanoma)
plot(fit2.aj,confint=FALSE,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
plot(fit2.arr,add=TRUE,col=1:5,lty=3,newdata=data.frame(logthick=quantile(Melanoma$logthick)))

## logistic risk regression
fit2.lrr <- LRR(Hist(time,status)~logthick,data=Melanoma,cause=1)
summary(fit2.lrr)

## change model for censoring weights
fit2a.lrr <- LRR(Hist(time,status)~logthick,
                 data=Melanoma,
                 cause=1,
                 cens.model="cox",
                 cens.formula=~sex+epicel+ulcer+age+logthick)
summary(fit2a.lrr)

##  compare prediction errors
\dontrun{
    library(pec)
    plot(pec(list(ARR=fit2.arr,AJ=fit2.aj,LRR=fit2.lrr),data=Melanoma,maxtime=3000))
}
# }}}

# {{{ multiple regression
library(pec)
library(riskRegression)
library(prodlim)
# absolute risk model
multi.arr <- ARR(Hist(time,status)~logthick+sex+age+ulcer,data=Melanoma,cause=1)

# stratified model allowing different baseline risk for the two gender
multi.arr <- ARR(Hist(time,status)~thick+strata(sex)+age+ulcer,data=Melanoma,cause=1)

# stratify by a continuous variable
multi.arr <- ARR(Hist(time,status)~tp(thick,power=1)+strata(sex)+age+ulcer,
                 data=Melanoma,
                 cause=1)

fit.arr2a <- ARR(Hist(time,status)~const(thick,power=1),data=Melanoma,cause=1)
fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
summary(fit.arr)

## logistic risk model
fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
summary(fit.lrr)

# }}}



## nearest neighbor non-parametric Aalen-Johansen estimate
library(prodlim)
fit.aj <- prodlim(Hist(time,status)~thick,data=Melanoma)
plot(fit.aj,confint=FALSE)

\dontrun{
    # prediction error
    library(pec)
    x <- pec(list(fit.arr2,fit.arr2a,fit.arr2b,fit.lrr),
             data=Melanoma,
             formula=Hist(time,status)~1,
             cause=1,
             B=10,
             splitMethod="none")
}
#'
}
\author{
Thomas H. Scheike \email{ts@biostat.ku.dk}

Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\references{
Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute
risk regression for competing risks: interpretation, link
functions and prediction Research report 11/8. Department
of Biostatistics, University of Copenhagen

Scheike, Zhang and Gerds (2008), Predicting
cumulativeincidence probability by direct binomial
regression,Biometrika, 95, 205-220.

Scheike and Zhang (2007), Flexible competing risks
regression modelling and goodness of fit, LIDA, 14,
464-483.

Martinussen and Scheike (2006), Dynamic regression models
for survival data, Springer.
}
\keyword{survival}

