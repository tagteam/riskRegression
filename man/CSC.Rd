\name{CSC}
\alias{CSC}
\title{
  Cause-specific Cox proportional hazard regression 
}
\description{
  Interface for fitting cause-specific Cox proportional hazard regression models in competing risk.
}
\usage{
CSC(formula, data, cause, survtype = "hazard", ...)
}
\arguments{
  \item{formula}{
    A list of formulae, one for each cause, each specifying a cause-specific
    Cox regression model. 
  }
  \item{data}{
    A data in which to fit the models.
  }
  \item{cause}{
    The cause of interest. Defaults to the first cause.
  }
  \item{survtype}{
    Either \code{"hazard"} (the default) or \code{"survival"}. If
    \code{"hazard"} fit cause-specific Cox regression models for all causes. If
    \code{"survival"} fit one cause-specific Cox regression model for the cause of
    interest and also a Cox regression model for event-free survival.
  }
  \item{\dots}{
   Arguments given to \code{coxph}. 
  }
}
% \details{
  
% }
\value{
%%  ~Describe the value returned
%  If it is a LIST, use
\item{call }{the call}
\item{models }{a list with the fitted (cause-specific) Cox regression
  objects}
\item{response }{the event history response }
\item{eventTimes }{the sorted (unique) event times }
\item{survtype }{the value of \code{survtype}}
\item{theCause }{the cause of interest. see \code{cause}}
\item{causes }{the other causes}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  Thomas A. Gerds
  \email{tag@biostat.ku.dk}
  
  Ulla B. Mogensen
  \email{ulmo@biostat.ku.dk}
}
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{coxph}}
}
\examples{
data(Melanoma)
## fit two cause-specific Cox models
## different formula for the two causes
fit1 <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),data=Melanoma)
print(fit1)

## same formula for both causes
fit2 <- CSC(Hist(time,status)~invasion+epicel+age,data=Melanoma)
print(fit2)

## combine a cause-specific Cox regression model for cause 2 
## and a Cox regression model for the event-free survival:
## different formula for cause 2 and event-free survival
fit3 <- CSC(list(Hist(time,status)~sex+invasion+epicel+age,Hist(time,status)~invasion+epicel+age),data=Melanoma)
print(fit3)

## same formula for both causes
fit4 <- CSC(Hist(time,status)~invasion+epicel+age,data=Melanoma,survtype="surv")
print(fit4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
