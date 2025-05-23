% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSC.R
\name{CSC}
\alias{CSC}
\title{Cause-specific Cox proportional hazard regression}
\usage{
CSC(formula, data, cause, surv.type = "hazard", fitter = "coxph", ...)
}
\arguments{
\item{formula}{Either a single \code{Hist} formula or a list of formulas.
If it is a list it must contain as many \code{Hist} formulas as there are
causes when \code{surv.type="hazard"} and exactly two formulas when \code{surv.type="survival"}.
If it is a list the first formula is used for the cause of interest specific Cox regression
and the other formula(s) either for the other cause specific Cox regression(s) or for the
Cox regression of the combined event where each cause counts as event. Note that when only one
formula is given the covariates enter in exactly the same way into all Cox regression analyses.}

\item{data}{A data in which to fit the models.}

\item{cause}{The cause of interest. Defaults to the first cause (see Details).}

\item{surv.type}{Either \code{"hazard"} (the default) or
\code{"survival"}.  If \code{"hazard"} fit cause-specific Cox
regression models for all causes.  If \code{"survival"} fit one
cause-specific Cox regression model for the cause of interest and
also a Cox regression model for event-free survival.}

\item{fitter}{Character string specifying the routine to fit the Cox regression models. Available are
\code{"coxph"} for \link[survival]{coxph}, \code{"coxph"} for \link[rms]{cph},
\code{"phreg"} for \link[mets]{phreg}, and \code{"glmnet"} for \link[glmnet]{glmnet}.}

\item{...}{Arguments given to the function defined by argument \code{fitter}.}
}
\value{
\item{models }{a list with the fitted (cause-specific) Cox
regression objects} \item{response }{the event history response }
\item{eventTimes }{the sorted (unique) event times } \item{surv.type }{the
value of \code{surv.type}} \item{theCause }{the cause of interest. see
\code{cause}} \item{causes }{the other causes} %% ...
}
\description{
Interface for fitting cause-specific Cox proportional hazard regression
models in competing risk.
}
\details{
The causes and their order
are determined by \code{prodlim::getStates()} applied to the Hist object.
}
\examples{

library(prodlim)
library(survival)
data(Melanoma)
## fit two cause-specific Cox models
## different formula for the two causes
fit1 <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
            data=Melanoma)
print(fit1)

\dontrun{
library(rms)
fit1a <- CSC(list(Hist(time,status)~sex+rcs(age,3),Hist(time,status)~invasion+epicel+log(thick)),
            data=Melanoma,fitter="cph")
print(fit1a)
}
\dontrun{
library(glmnet)
# lasso
fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
            data=Melanoma,fitter="glmnet")
# rigde regression
fit1c <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
            data=Melanoma,fitter="glmnet")
print(fit1b)
}
\dontrun{
library(Publish)
publish(fit1)
}

## model hazard of all cause mortality instead of hazard of type 2
fit1a <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
             data=Melanoma,
             surv.type="surv")

## the predicted probabilities are similar
plot(predictRisk(fit1,times=500,cause=1,newdata=Melanoma),
     predictRisk(fit1a,times=500,cause=1,newdata=Melanoma))

## special case where cause 2 has no covariates
fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~1),
             data=Melanoma)
print(fit1b)
predict(fit1b,cause=1,times=100,newdata=Melanoma)


## same formula for both causes
fit2 <- CSC(Hist(time,status)~invasion+epicel+age,
            data=Melanoma)
print(fit2)

## combine a cause-specific Cox regression model for cause 2
## and a Cox regression model for the event-free survival:
## different formula for cause 2 and event-free survival
fit3 <- CSC(list(Hist(time,status)~sex+invasion+epicel+age,
                 Hist(time,status)~invasion+epicel+age),
            surv.type="surv",
            data=Melanoma)
print(fit3)

## same formula for both causes
fit4 <- CSC(Hist(time,status)~invasion+epicel+age,
            data=Melanoma,
            surv.type="surv")
print(fit4)

## strata
fit5 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
            data=Melanoma,
            surv.type="surv")
print(fit5)

## sanity checks

cox1 <- coxph(Surv(time,status==1)~invasion+epicel+age+strata(sex),data=Melanoma)
cox2 <- coxph(Surv(time,status!=0)~invasion+epicel+age+strata(sex),data=Melanoma)
all.equal(coef(cox1),coef(fit5$models[[1]]))
all.equal(coef(cox2),coef(fit5$models[[2]]))

## predictions
##
## surv.type = "hazard": predictions for both causes can be extracted
## from the same fit
fit2 <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma)
predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predictRisk(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predictRisk(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
predict(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))

## surv.type = "surv" we need to change the cause of interest
library(survival)
fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
            data=Melanoma,
            surv.type="surv",cause=2)
## now this does not work because the object was fitted with surv.type='surv'
try(predictRisk(fit5.2,cause=1,newdata=Melanoma,times=4000))

## but this does
predictRisk(fit5.2,cause=2,newdata=Melanoma,times=100)
predict(fit5.2,cause=2,newdata=Melanoma,times=100)
predict(fit5.2,cause=2,newdata=Melanoma[4,],times=100)

fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
            data=Melanoma,
            surv.type="hazard",cause=2)

}
\references{
B. Ozenne, A. L. Soerensen, T.H. Scheike, C.T. Torp-Pedersen,
and T.A. Gerds. riskregression: Predicting the risk
of an event using Cox regression models. R Journal, 9(2):440--460, 2017.

J Benichou and Mitchell H Gail. Estimates of absolute cause-specific risk
in cohort studies. Biometrics, pages 813--826, 1990.

T.A. Gerds, T.H. Scheike, and P.K. Andersen. Absolute risk regression for
competing risks: Interpretation, link functions, and prediction. Statistics
in Medicine, 31(29):3921--3930, 2012.
}
\seealso{
\code{\link[survival]{coxph}}
}
\author{
Thomas A. Gerds \email{tag@biostat.ku.dk} and Ulla B. Mogensen
}
\keyword{survival}
