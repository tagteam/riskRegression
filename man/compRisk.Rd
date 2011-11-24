\name{riskRegression}
\alias{riskRegression}
\title{Competings Risks Regression} 
\description{
Fits a semiparametric regression model for the cumulative incidence function:
\deqn{
P(T < t, cause=1 | x,z) = P_1(t,x,z) = h(g(t,x,z))
}
for a known link-function \eqn{h()},
predictor variables with time-constant effect z,
predictor variables with time-varying effect x, 
and known prediction-function \eqn{g(t,x,z)}.
for the probability of dying from cause 1 in a situation with 
competing causes of death. 

We consider the following models :

(default) the "relative" model where
\eqn{h(x)=\exp(-x)} and 
\deqn{
  g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}

The interpretation of the regression coefficients 

(optional) the additive model where
\eqn{h(x)=1-\exp(-x)} and 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta

}
2)  
the proportional setting that includes the Fine & Gray (FG) model and
some extensions where 
\eqn{h(x)=1-\exp(-\exp(x))} and 
\deqn{
g(t,x,z) = \exp(x^T A(t) + (diag(t^p) z)^T \beta)
}
The FG model is obtained when \eqn{x=1}. 

3) a logistic model where 
\eqn{h(x)=\exp(x)/( 1+\exp(x))} and 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}

Where p by default is 1 for the additive model and 0 for 
the other models.
In general p may be powers of the same length as z. 
}
\usage{
riskRegression(formula,
	 data = sys.parent(),
	 times = NULL,
	 link = "relative",
	 cause = 1,
	 resample.iid = 1,
	 cens.model = "KM",
	 n.sim = 500,
	 maxiter = 50,
	 silent = 1,
	 convLevel = 6,
	 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    Example: Hist(time,status)~age + const(sex,power=1)
  }
  \item{data}{
  }
  \item{times}{
    Vector of times. For each time point in \code{times} estimate the baseline risk and the timevarying
    coefficients.
  }
  \item{link}{
  }
  \item{cause}{
    The cause of interest. 
  }
  \item{resample.iid}{
%%     ~~Describe \code{resample.iid} here~~
}
  \item{cens.model}{
%%     ~~Describe \code{cens.model} here~~
}
  \item{n.sim}{
%%     ~~Describe \code{n.sim} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{silent}{
%%     ~~Describe \code{silent} here~~
}
  \item{convLevel}{
%%     ~~Describe \code{convLevel} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
