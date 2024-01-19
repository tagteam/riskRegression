### test-predictRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 16 2022 (08:10) 
## Version: 
## Last-Updated: Dec 19 2023 (08:34) 
##           By: Thomas Alexander Gerds
##     Update #: 7
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

## OBS: Changed the class so that they are not two definitions of predictRisk.
GrpSurv2 <- function(formula,data,...){
    requireNamespace("grpreg")
    EHF = EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
    fit = grpsurv(X = EHF$design,y = EHF$event.history,...)
    fit = list(grpsurv.object = fit,terms = terms(formula),call=match.call())
    class(fit) = c("GrpSurv2",class(fit))
    fit
}

predictRisk.GrpSurv2 <- function(object, newdata, times, cause, ...){
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
f1 <- GrpSurv2(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,lambda = .1,group = c(1,1,rep(2,8)))
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
f1 <- Hal9001(Surv(time,event)~X1+X6+X7+X9,data =d)
f2 <- coxph(Surv(time,event)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data = d,x = 1,y = 1)
p <- predictRisk(object = f1,newdata = nd,times = c(5,10))
x <- Score(list(hal = f1,cox = f2),
           formula = Surv(time,event)~1,
           data = nd,
           times = c(5,10))
summary(x)

library(riskRegression)
library(survival)

set.seed(245)

data(nafld, package="survival")

# generating random cluster
nafld1$clust <- sample(1:10, size=nrow(nafld1), replace=TRUE)
nafld1$cluster <- nafld1$clust

# fitting a model with a frailty term without naming it "cluster"
model1 <- coxph(Surv(futime, status) ~ age + male + frailty(clust),
               data=nafld1, x=TRUE)

# fitting a model with a frailty term named "cluster"
model2 <- coxph(Surv(futime, status) ~ age + male + frailty(cluster),
                data=nafld1, x=TRUE)


# getting predictions
test1 <- predictRisk(model1, nafld1, times=100)
test2 <- predictRisk(model2, nafld1, times=100) 

all.equal(test1, test2)

######################################################################
### test-predictRisk.R ends here
