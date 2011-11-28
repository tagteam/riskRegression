\name{Melanoma}
\alias{Melanoma}
\docType{data}
\title{
  Malignant melanoma data
}
\description{
In the period 1962-77, 205 patients with malignant melanoma (cancer of
the skin) had a radical operation performed at Odense University
Hospital, Denmark. All patients were followed until the end of 1977 by 
which time 134 were still alive while 71 had died (of out whom 57 had
died from cancer and 14 from other causes).

The object of the study was to assess the effect of risk factors on
survival. Among such risk factors were the sex and age of the patients
and the histological variables tumor thickness and ulceration (absent
vs. present).

}
\usage{data(Melanoma)}
\format{
  A data frame with 205 observations on the following 12 variables.
  \describe{
    \item{\code{time}}{ time in days from operation}
    \item{\code{event}}{a factor with levels \code{censored}
      \code{death.malignant.melanoma} \code{death.other.causes}}
    \item{\code{status}}{a numeric with values \code{0=censored} \code{1=death.malignant.melanoma} \code{2=death.other.causes}}
    \item{\code{invasion}}{a factor with levels \code{level.0}, \code{level.1}, \code{level.2}}
    \item{\code{ici}}{inflammatory cell infiltration (ICI): 0, 1, 2 or 3}
    \item{\code{epicel}}{a factor with levels \code{not present} \code{present}}
    \item{\code{ulcer}}{a factor with levels \code{not present} \code{present}}
    \item{\code{thick}}{tumour thickness (in 1/100 mm)}
    \item{\code{sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{age}}{age at operation (years)}
  }
}
\source{
\url{http://192.38.117.59/~linearpredictors/?page=datasets&dataset=Melanoma}
}
\references{
  Regression with linear predictors (2010)
  
  Andersen, P.K. and Skovgaard, L.T.
  
  Springer Verlag
}
\examples{
data(Melanoma)
}
\keyword{datasets}
