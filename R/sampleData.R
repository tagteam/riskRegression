### sampleData.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (09:43) 
## Version: 
## last-updated: Sep  4 2017 (11:12) 
##           By: Thomas Alexander Gerds
##     Update #: 32
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Simulate data with binary outcome and 10 covariates.
##'
##' For the actual lava::regression parameters see the function definition.
##' @title Simulate data with binary or time-to-event outcome
##' @aliases sampleDataTD
##' @param n Sample size
##' @param n.intervals \code{sampleDataTD} only: the maximum number of episodes in which the covariates are updated.
##' @param outcome Character vector. Response variables are generated
##' according to keywords: \code{"binary"} = binary response,
##' \code{"survival"} = survival response, \code{"competing.risks"} =
##' competing risks response
##' @param formula Specify regression coefficients
##' @return Simulated data as data.table with n rows and the following columns:
##' Y (binary outcome), time (non-binary outcome), event (non-binary outcome), X1-X5 (binary predictors), X6-X10 (continous predictors)
##' @seealso lvm
##' @examples
##' sampleData(10,outcome="binary")
##' sampleData(10,outcome="survival")
##' sampleData(10,outcome="competing.risks")
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
sampleData <- function(n,outcome="competing.risks",formula= ~ f(X1,2) + f(X2,-0.033) + f(X3,0.4) + f(X6,.1) + f(X7,-.1) + f(X8,.5) + f(X9,-1)){
    outcome <- match.arg(outcome,c("survival","competing.risks","binary"))
    m <- lava::lvm()
    lava::distribution(m,~X6) <- lava::normal.lvm(mean=60,sd=15)
    lava::distribution(m,~X7) <- lava::normal.lvm(mean=60,sd=5)
    lava::distribution(m,~X8) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X9) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X10) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X1) <- lava::binomial.lvm(p=c(.1))
    lava::distribution(m,~X2) <- lava::binomial.lvm(p=c(.2))
    lava::distribution(m,~X3) <- lava::binomial.lvm(p=c(.3))
    lava::distribution(m,~X4) <- lava::binomial.lvm(p=c(.4))
    lava::distribution(m,~X5) <- lava::binomial.lvm(p=c(.5))
    if ("binary"%in%outcome){
        lava::distribution(m,~Y) <- lava::binomial.lvm()
        lava::regression(m) <- stats::update(formula,"Y~.")
    }
    if ("survival"%in%outcome){
        lava::distribution(m, "eventtime") <- lava::coxWeibull.lvm(scale = 1/100)
        lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
        m <- lava::eventTime(m, time ~ min(eventtime = 1, censtime = 0),"event")
        lava::regression(m) <- stats::update(formula,"eventtime~.")
    }
    if ("competing.risks"%in%outcome){
        lava::distribution(m, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/100)
        lava::distribution(m, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/100)
        lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
        m <- lava::eventTime(m, time ~ min(eventtime1 = 1, eventtime2 = 2, censtime = 0), "event")
        lava::regression(m) <- stats::update(formula,"eventtime1~.")
    }
    data.table::as.data.table(lava::sim(m,n))
}

##' @export 
sampleDataTD <- function(n,n.intervals=5,outcome="competing.risks",formula= ~ f(X1,2) + f(X2,-0.033) + f(X3,0.4) + f(X6,.1) + f(X7,-.1) + f(X8,.5) + f(X9,-1)){
    m <- lava::lvm()
    lava::distribution(m,~X6) <- lava::normal.lvm(mean=60,sd=15)
    lava::distribution(m,~X7) <- lava::normal.lvm(mean=60,sd=5)
    lava::distribution(m,~X8) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X9) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X10) <- lava::normal.lvm(mean=0,sd=1)
    lava::distribution(m,~X1) <- lava::binomial.lvm(p=c(.1))
    lava::distribution(m,~X2) <- lava::binomial.lvm(p=c(.2))
    lava::distribution(m,~X3) <- lava::binomial.lvm(p=c(.3))
    lava::distribution(m,~X4) <- lava::binomial.lvm(p=c(.4))
    lava::distribution(m,~X5) <- lava::binomial.lvm(p=c(.5))
    lava::distribution(m, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/100)
    lava::distribution(m, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/100)
    lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/10)
    m <- lava::eventTime(m, time ~ min(eventtime1 = 1, eventtime2 = 2, censtime = 0), "event")
    lava::regression(m) <- stats::update(formula,"eventtime1~.")
    update.mydata <- function(m,data){
        is.censored <- data$event==0
        d <- lava::sim(m,sum(is.censored))
        data.table::setDT(d)
        d[,start:=data[is.censored,time]]
        if (any(is.na(d$start))) browser()
        d[,time:=time+data[is.censored,time]]
        d
    }
    dTD <- vector(n.intervals,mode="list")
    dTD[[1]] <- data.table::as.data.table(lava::sim(m,n))
    dTD[[1]][,start:=0]
    for (t in 2:n.intervals){
        dTD[[t]] <- update.mydata(m,dTD[[t-1]])
    }
    rbindlist(dTD)
}

#----------------------------------------------------------------------
### sampleData.R ends here
