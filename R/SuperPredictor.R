### SuperPredictor.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  6 2019 (18:22) 
## Version: 
## Last-Updated: Mar  9 2022 (15:49) 
##           By: Thomas Alexander Gerds
##     Update #: 10
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Formula interface for SuperLearner::SuperLearner
##'
##' Formula interface for SuperLearner::SuperLearner
##' ##' @param formula
##' @title Formula interface for SuperLearner::SuperLearner
##' @param formula where the left hand side specifies the outcome and the right hand side the predictors
##' @param data data set in which formula can be evaluated
##' @param family the outcome family. default is binomial
##' @param SL.library the SuperLearner libraries
##' @param ... passed to SuperLearner::SuperLearner
##' @examples
##' \dontrun{
##' if(require("SuperLearner",quietly=TRUE)){
##' library(SuperLearner)
##' library(data.table)
##' set.seed(10)
##' d = sampleData(338, outcome="binary")
##' spfit = SuperPredictor(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=d)
##' predictRisk(spfit)
##' x <- Score(list(spfit),data=d,formula=Y~1)
##' }
##' }
##' @export 
SuperPredictor <- function(formula,data,family="binomial",SL.library=c("SL.glm","SL.glm.interaction","SL.ranger"),...){
    vv <- all.vars(formula)
    yy <- vv[[1]]
    xx <- vv[-1]
    if (is.data.table(data)) XX <- data[,xx,with=FALSE] else XX <- data[,xx]
    YY <- data[[yy]]
    if (is.factor(YY)) YY <- as.numeric(YY!=levels(YY)[1])
    fit <- SuperLearner::SuperLearner(Y=YY,
                                      X=XX,
                                      family=family,
                                      SL.library=SL.library,
                                      ...)
    class(fit) <- "SuperPredictor"
    fit$call <- match.call()
    fit
}

######################################################################
### SuperPredictor.R ends here
