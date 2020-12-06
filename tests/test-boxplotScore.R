### test-boxplotScore.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Aug 23 2016 (17:07) 
## Version: 
## last-updated: Dec  5 2020 (12:38) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(testthat)
library(riskRegression)
library(data.table)
context("Retrospective boxplots of predicted risks")
# {{{ "boxplot.Score"
test_that("boxplot.Score",{
    d <- data.table(time=1:10,
                    Y=c(rep(1,5),rep(0,5)),
                    status=1,
                    b=c(.3,.2,.2,.3,.4),
                    a=c(.8,.5,.2,.2,0))
    d[,diffba:=b-a]
    ## d <- rbindlist(list(d,d,d))
    Qe <- d[time<5.1,quantile(diffba,probs=c(0.05,0.25,0.5,0.75,0.95),type=1)]
    Qef <- d[time>=5.1,quantile(diffba,probs=c(0.05,0.25,0.5,0.75,0.95),type=1)]

    ## [ERROR: * cannot open the connection]
    binCase <- Score(list(a=d$a,b=d$b),formula=Y~1,data=d,summary="riskQuantile",null.model=FALSE,plots=NULL,metrics=NULL)
    survCase <- Score(list(a=d$a,b=d$b),formula=Surv(time,status)~1,times=5.1,data=d,summary="riskQuantile",probs=c(0.05,0.25,0.5,0.75,0.95),null.model=FALSE,plots=NULL,metrics=NULL)
    expect_equal(as.numeric(unlist(survCase$riskQuantile$contrasts[cause=="event",5:9,with=FALSE])),as.numeric(Qe))
    expect_equal(as.numeric(unlist(survCase$riskQuantile$contrasts[cause=="event-free",5:9,with=FALSE])),as.numeric(Qef))    
})

# }}}
# {{{ "getQuantile"
test_that("getQuantile",{
    x <- 0:5
    Fx <- cumsum(rep(1,length(x)))/length(x)
    expect_equal(riskRegression:::getQuantile(x=x,Fx=Fx,Q=seq(0,1,.25)),
                 as.numeric(quantile(0:5,type=1,probs=seq(0,1,.25))))
    qseq <- runif(8,0,1)
    expect_equal(riskRegression:::getQuantile(x=x,Fx=Fx,Q=qseq),
                 as.numeric(quantile(0:5,type=1,probs=qseq)))
})
# }}}
#----------------------------------------------------------------------
### test-boxplotScore.R ends here
