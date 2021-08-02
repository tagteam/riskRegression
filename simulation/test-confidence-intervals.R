### test-confidence-intervals.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 20 2021 (16:19) 
## Version: 
## Last-Updated: Aug  2 2021 (13:49) 
##           By: Thomas Alexander Gerds
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary:
##
##  Working document for the collaboration between
##  Johan Sebastian Ohlendorff and Thomas Alexander Gerds.
##
##  In this document we evaluate the small sample properties of our
##  prediction performance estimators and their standard errors in
##  synthesized data. 
##
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(lava)
library(survival)
library(riskRegression)
data(pbc)
pbc <- na.omit(pbc)
pbc$treat <- 1*(pbc$trt==2)
    
# synthesize binary outcome data
mb <- synthesize(treat~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

# synthesize survival outcome data
ms <- synthesize(Hist(time,status!=0)~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

# synthesize competing risk outcome data
mc <- synthesize(Hist(time,status)~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

## Setting 1: learn/test
## --------------------------------------------------------------------

### A. binary  outcome

## We have a dataset for building the risk prediction model (learndata)
## and a second dataset for estimating the prediction performance (testdata).
# 
## We want to study the small-sample behaviour of the standard errors as obtained with
## the Score function and the derived confidence intervals and p-values.
# 
## In this setting the learndata are fixed
## and the simulation runs across B test datasets.
#
set.seed(17)
largedata <- sim(mb,30000)
set.seed(8)
learndata <- sim(mb,300)
model1 <- glm(treat~logbili+logprotime+sex+age,data=learndata,family="binomial")
model2 <- glm(treat~logprotime+sex+age,data=learndata,family="binomial")
# find true value (may need to be repeated)
true.X <- Score(list(model1,model2),data=largedata,formula=treat~1,metrics="brier")
# run one
testdata <- sim(mb,300)
x <- Score(list(m1=model1,m2=model2),data=testdata,formula=treat~1,metrics="brier",null.model=0L)
x$Brier$score
x$Brier$contrast

# run simulation
run <- function(...,n=300){
    testdata <- sim(mb,n=n)
    x <- Score(list(m1=model1,m2=model2),
               data=testdata,
               formula=treat~1,
               metrics="brier",null.model=0L)
    # extract standard error and confidence intervals for the 3 models
    c(Brier=x$Brier$score$Brier,
      se.Brier=x$Brier$score$se,
      lower.Brier=x$Brier$score$lower,
      upper.Brier=x$Brier$score$upper,
      # comparison of models
      delta.Brier=x$Brier$contrast$delta.Brier,
      se.delta.Brier=x$Brier$contrast$se,
      lower.delta.Brier=x$Brier$contrast$lower,
      upper.delta.Brier=x$Brier$contrast$upper)
}

## repeat ones
run(1)
# same as
sim(run,1)

# repeat 100 times (should be 10,000 times)
sim.result <- sim(run,100)

# get approximate true values
true.values <- c(Brier=true.X$Brier$score$Brier,
                 se.Brier=true.X$Brier$score$se,
                 lower.Brier=true.X$Brier$score$lower,
                 upper.Brier=true.X$Brier$score$upper,
                 # comparison of models
                 delta.Brier=true.X$Brier$contrast$delta.Brier,
                 se.delta.Brier=true.X$Brier$contrast$se,
                 lower.delta.Brier=true.X$Brier$contrast$lower,
                 upper.delta.Brier=true.X$Brier$contrast$upper)


summary(sim.result,true=true.values)

### B. survival outcome

### C. competing risk outcome



## Setting 2: leave-one-out bootstrap
## --------------------------------------------------------------------

## Setting 3: cross-validation
## --------------------------------------------------------------------


#----------------------------------------------------------------------
### test-confidence-intervals.R ends here
