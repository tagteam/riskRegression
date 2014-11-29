\name{predict.riskRegression}
\alias{predict.riskRegression}
\title{Predict individual risk.}
\usage{
\method{predict}{riskRegression}(object, newdata, ...)
}
\arguments{
  \item{object}{Fitted object obtained with one of
  \code{ARR}, \code{LRR}, \code{riskRegression}.}

  \item{newdata}{A data frame containing predictor variable
  combinations for which to compute predicted risk.}

  \item{\dots}{not used}
}
\description{
Extract predictions from a risk prediction model.
}
\examples{
data(Melanoma)
library(prodlim)
fit.tarr <- ARR(Hist(time,status)~age+invasion+strata(sex),data=Melanoma,cause=1)
predict(fit.tarr,newdata=data.frame(age=48,
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Female",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=data.frame(age=48,
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Male",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=data.frame(age=c(48,58,68),
                     invasion=factor("level.1",
                         levels=levels(Melanoma$invasion)),
                     sex=factor("Male",levels=levels(Melanoma$sex))))
predict(fit.tarr,newdata=Melanoma[1:4,])
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
}
\keyword{survival}

