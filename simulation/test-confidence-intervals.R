### test-confidence-intervals.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 20 2021 (16:19) 
## Version: 
## Last-Updated: Jul 23 2021 (10:33) 
##           By: Thomas Alexander Gerds
##     Update #: 6
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

# synthesize binary outcome data
mb <- synthesize(treat~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

# synthesize survival outcome data
ms <- synthesize(Hist(time,status!=0)~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

# synthesize competing risk outcome data
mc <- synthesize(Hist(time,status)~log(bili)+log(protime)+edema+sex+age,data=pbc,na.rm=TRUE)

## Setting 1: learn/test
## --------------------------------------------------------------------

## We have a dataset for building the risk prediction model (learndata)
## and a second dataset for estimating the prediction performance (testdata).
learndata <- 
x <- Score






#----------------------------------------------------------------------
### test-confidence-intervals.R ends here
