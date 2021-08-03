### develop-synthesizer.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 21 2021 (08:59)
## Version:
## Last-Updated: Aug  2 2021 (13:26)
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
##
### Commentary:
##
##  Working document for the collaboration between
##  Johan Sebastian Ohlendorff and Thomas Alexander Gerds.
##
## In this document we develop the functionality of the data synthesizer.
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

library(riskRegression)
library(survival)
library(lava)
data(pbc)
pbc <- na.omit(pbc)

# this works
s <- synthesize(Surv(time,status)~age+sex+protime,data=pbc,verbose=FALSE)
sim(s,10)

# Task 1: learning categorical variables
# remaining problem: line 310 in synthesize.R
#  reg_formula = time.event.0 ~ age + sex + edema + protime
#  coef(G)[-1]/G$scale = age sexf edema0.5 edema1 protime
#  there are 2 coefficients, one for each contrast defined by the categorical variable
#  edema=0.5/edema=0 and one for edema=1/edema=0.
#  there are several ideas to fix this without the dummy
#  variables (p-1 dummy 0-1 variables are used to parametrize a categorical with p levels).
#  the dummy variables are created in lines 238-239 of synthesize.R
#  to fix the problem quickly, we can try to modify object$M such that the
#  regressions of the time variables are conditioning on edema0.5 and
#  edema1 instead of edema.
pbc$event <- 1*(pbc$status!=0)
synthesize(Surv(time,event)~sex+edema,data=pbc)

s <- synthesize(Surv(time,status)~age+sex+edema+protime,data=pbc)
sim(s,10)

# Task 2: learning log transformations
# this should work
s <- synthesize(Surv(time,event)~age+sex+log(protime),data=pbc)
d <- sim(s,1000)
d$time <- as.numeric(d$time)
d$event <- 1*(d$status!=0)

fit <- coxph(Surv(time,event)~age+sex+log(protime),data=pbc)
fit.s <- coxph(Surv(time,event)~age+sex+logprotime,data=d)
cbind(coef(fit),coef(fit.s))

# Task 3: learning structural equations from data
# when (new option) recursiv is TRUE we want that the following formula is
# read such that we have regression equations that model the dependence
# between the covariates:
## age~sex
## edema~age+sex
## protime~edema+age+sex
s <- synthesize(Surv(time,status)~protime+age+sex,data=pbc,recursive=TRUE)

d <- sim(s,1000)
d$time <- as.numeric(d$time)
d$event <- 1*(d$status!=0)

fit.protime <- glm(protime~age+sex+protime,data=pbc)
fit.protime.s <- glm(protime~age+sex+protime,data=d)
cbind(coef(fit.protime),coef(fit.protime.s))




#----------------------------------------------------------------------
### develop-synthesizer.R ends here
