### develop-synthesizer.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 21 2021 (08:59) 
## Version: 
## Last-Updated: Jul 21 2021 (10:16) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
s <- synthesize(Surv(time)~age+sex+protime,data=pbc,verbose=FALSE)
sim(s,10)

# Task 1: learning categorical variables
# this should work
s <- synthesize(Surv(time)~age+sex+edema+protime,data=pbc)

# Task 2: learning log transformations
# this should work 
s <- synthesize(Surv(time)~age+sex+edema+log(protime),data=pbc)

# Task 3: learning structural equations from data
# when (new option) recursiv is TRUE we want that the following formula is
# read such that we have regression equations that model the dependence
# between the covariates:
## age~sex
## edema~age+sex
## protime~edema+age+sex
s <- synthesize(Surv(time)~protime+edema+age+sex,data=pbc,recursive=TRUE)



#----------------------------------------------------------------------
### develop-synthesizer.R ends here
