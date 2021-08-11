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
pbc$event <- 1*(pbc$status!=0)

# synthesize binary outcome data
mb <- synthesize(treat~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

# synthesize survival outcome data
ms <- synthesize(Hist(time,event)~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

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
Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2])
Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])

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
      upper.delta.Brier=x$Brier$contrast$upper,
      # coverage
      Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2]),
      Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])
      )
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
# WARNING!!! This does not work! R crashes!!! As you have correctly specified!
set.seed(17)
largedata <- sim(ms,30000)
set.seed(8)
learndata <- sim(ms,300)
model1 <- coxph(Surv(time,event)~logbili+logprotime+sex+age,data=learndata)
model2 <- coxph(Surv(time,event)~logprotime+sex+age,data=learndata)
# find true value (may need to be repeated)
# correct formula?
true.X <- Score(list(model1,model2),data=largedata,formula=Surv(time,event)~1,metrics="brier")

# run simulation
run <- function(...,n=300){
  testdata <- sim(mb,n=n)
  x <- Score(list(m1=model1,m2=model2),
             data=testdata,
             formula=Surv(time,event)~1,
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
    upper.delta.Brier=x$Brier$contrast$upper,
    # coverage
    Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2]),
    Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])
  )
}

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

### C. competing risk outcome
# WARNING!!! This does not work! R crashes!!! As you have correctly specified!
set.seed(17)
largedata <- sim(mc,30000)
set.seed(8)
learndata <- sim(mc,300)
model1 <- coxph(Surv(time,status)~logbili+logprotime+sex+age,data=learndata)
model2 <- coxph(Surv(time,status)~logprotime+sex+age,data=learndata)
# find true value (may need to be repeated)
# correct formula?
true.X <- Score(list(model1,model2),data=largedata,formula=Surv(time,status)~1,metrics="brier")

# run simulation
run <- function(...,n=300){
  testdata <- sim(mb,n=n)
  x <- Score(list(m1=model1,m2=model2),
             data=testdata,
             formula=Surv(time,event)~1,
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
    upper.delta.Brier=x$Brier$contrast$upper,
    # coverage
    Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2]),
    Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])
  )
}

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


## Setting 2: leave-one-out bootstrap
## --------------------------------------------------------------------
# done for the logistic regression, but can probably be similarly done for survival
set.seed(17)
largedata <- sim(mb,30000)
set.seed(8)
learndata <- sim(mb,300)
model1 <- glm(treat~logbili+logprotime+sex+age,data=learndata,family="binomial")
model2 <- glm(treat~logprotime+sex+age,data=learndata,family="binomial")
# find true value (may need to be repeated)
true.X <- Score(list(model1,model2),data=largedata,formula=treat~1,metrics="brier")

# run simulation
run <- function(...,n=300){
  testdata <- sim(mb,n=n)
  x <- Score(list(m1=model1,m2=model2),
             data=testdata,
             formula=treat~1,
             metrics="brier",split.method="loob",B=200,null.model=0L)
  # extract standard error and confidence intervals for the 3 models
  c(Brier=x$Brier$score$Brier,
    se.Brier=x$Brier$score$se,
    lower.Brier=x$Brier$score$lower,
    upper.Brier=x$Brier$score$upper,
    # comparison of models
    delta.Brier=x$Brier$contrast$delta.Brier,
    se.delta.Brier=x$Brier$contrast$se,
    lower.delta.Brier=x$Brier$contrast$lower,
    upper.delta.Brier=x$Brier$contrast$upper,
    # coverage
    Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2]),
    Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])
  )
}

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


## Setting 3: cross-validation
## --------------------------------------------------------------------
# done for the logistic regression, but can probably be similarly done for survival
set.seed(17)
largedata <- sim(mb,30000)
set.seed(8)
learndata <- sim(mb,300)
model1 <- glm(treat~logbili+logprotime+sex+age,data=learndata,family="binomial")
model2 <- glm(treat~logprotime+sex+age,data=learndata,family="binomial")
# find true value (may need to be repeated)
true.X <- Score(list(model1,model2),data=largedata,formula=treat~1,metrics="brier")

# run simulation
run <- function(...,n=300){
  testdata <- sim(mb,n=n)
  x <- Score(list(m1=model1,m2=model2),
             data=testdata,
             formula=treat~1,
             metrics="brier",split.method="bootcv",B=100,null.model=0L)
  # extract standard error and confidence intervals for the 3 models
  c(Brier=x$Brier$score$Brier,
    se.Brier=x$Brier$score$se,
    lower.Brier=x$Brier$score$lower,
    upper.Brier=x$Brier$score$upper,
    # comparison of models
    delta.Brier=x$Brier$contrast$delta.Brier,
    se.delta.Brier=x$Brier$contrast$se,
    lower.delta.Brier=x$Brier$contrast$lower,
    upper.delta.Brier=x$Brier$contrast$upper,
    # coverage
    Cov1=1*(x$Brier$score$lower[1]<=true.X$Brier$score$Brier[2] && x$Brier$score$upper[1]>=true.X$Brier$score$Brier[2]),
    Cov2=1*(x$Brier$score$lower[2]<=true.X$Brier$score$Brier[3] && x$Brier$score$upper[2]>=true.X$Brier$score$Brier[3])
  )
}

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

#----------------------------------------------------------------------
### test-confidence-intervals.R ends here
