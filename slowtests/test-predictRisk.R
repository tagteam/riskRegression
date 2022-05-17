### test-predictRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 16 2022 (08:10) 
## Version: 
## Last-Updated: May 17 2022 (13:52) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(grpreg)
library(survival)
library(prodlim)
library(riskRegression)
GrpSurv <- function(formula,data,...){
    requireNamespace("grpreg")
    EHF = EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
    fit = grpsurv(X = EHF$design,y = EHF$event.history,...)
    fit = list(grpsurv.object = fit,terms = terms(formula),call=match.call())
    class(fit) = c("GrpSurv",class(fit))
    fit
}

predictRisk.GrpSurv <- function(object, newdata, times, cause, ...){
    newdata$dummy.time=rep(1,NROW(newdata))
    newdata$dummy.event=rep(1,NROW(newdata))
    rhs <- as.formula(delete.response(object$terms))
    dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                       data=newdata,
                                       specials = NULL,
                                       unspecialsDesign=TRUE)
    newdata$dummy.time = NULL
    newdata$dummy.event = NULL
    p <- predict(object$grpsurv.object, EHF$design, type="survival")
    p <- 1-sapply(p,function(f)f(times))
    if (length(times) == 1)
        p = cbind(p)
    else(p = t(p))
    p
}


set.seed(8)
d <- sampleData(200,outcome = "survival")
nd <- sampleData(71,outcome = "survival")
f1 <- GrpSurv(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,lambda = .1,group = c(1,1,rep(2,8)))
f2 <- coxph(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,x = 1,y = 1)
p <- predictRisk(object = f1,newdata = nd,times = c(5,10))
x <- Score(list(grpsurv = f1,cox = f2),
           formula = Surv(time,event)~1,
           data = nd,
           times = c(5,10))
summary(x)


set.seed(8)
library(hal9001)
d <- sampleData(200,outcome = "survival")
nd <- sampleData(71,outcome = "survival")
f1 <- hal(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,lambda = .1,group = c(1,1,rep(2,8)))
f2 <- coxph(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,x = 1,y = 1)
p <- predictRisk(object = f1,newdata = nd,times = c(5,10))
x <- Score(list(grpsurv = f1,cox = f2),
           formula = Surv(time,event)~1,
           data = nd,
           times = c(5,10))
summary(x)



######################################################################
### test-predictRisk.R ends here
