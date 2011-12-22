\name{predict.riskRegression}
\alias{predict.riskRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict individual risk.
}
\description{
Extract predictions from a risk prediction model.
}
\usage{
\method{predict}{riskRegression}(object,
		       newdata,
		       ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Fitted object obtained with one of \code{ARR}, \code{LRR},
    \code{riskRegression}.
}
  \item{newdata}{
    A data frame containing predictor variable combinations for
   which to compute predicted risk.}
\item{\dots}{
  not used
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
 \references{
Gerds, TA and Scheike, T and Andersen, PK (2011)
Absolute risk regression for competing risks: interpretation, link functions and prediction
Research report 11/8. Department of Biostatistics, University of
Copenhagen
}
\author{
  Thomas H. Scheike
  \email{ts@biostat.ku.dk}
  
  Thomas A. Gerds
  \email{tag@biostat.ku.dk} 
}
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
data(Melanoma)

fit.tarr <- ARR(Hist(time,status)~age+invasion+strata(sex),data=Melanoma,cause=1)
predict(fit.tarr,newdata=data.frame(age=48,invasion="level.1",sex="Female"))
predict(fit.tarr,newdata=data.frame(age=48,invasion="level.1",sex="Male"))
predict(fit.tarr,newdata=data.frame(age=c(48,58,68),invasion="level.1",sex="Male"))
predict(fit.tarr,newdata=Melanoma[1:4,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival }
