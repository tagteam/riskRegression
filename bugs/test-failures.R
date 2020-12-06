### test-failures.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Dec  6 2020 (08:39) 
## Version: 
## Last-Updated: Dec  6 2020 (08:59) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(riskRegression)
library(testthat)
library(rms)
library(survival)
library(data.table)
library(lava)
test_that("ate in presence of splines",{
    ## ** Splines
    ## Tag: 10/06/20 6:21 
    n <- 500
    set.seed(10)
    dt <- sampleData(n, outcome = "competing.risks")
    dt$X1 <- factor(rbinom(n, prob = c(0.4) , size = 1), labels = paste0("T",0:1))
    m.event <-  CSC(Hist(time,event)~ X1+X2+X6+rms::rcs(X7,3),data=dt)
    ate.1 <- ate(event = m.event,
                 treatment = "X1",
                 data = dt, times = 1:5, 
                 cause = 1, product.limit = FALSE,
                 verbose = TRUE)
    suppressWarnings(ate.2 <- ate(event = m.event,
                                  treatment = "X1",
                                  data = dt, times = 1:5, 
                                  cause = 1, product.limit = FALSE,
                                  verbose = TRUE, B = 100))
    expect_output(print(ate.1))
    expect_output(print(ate.2))
})

## ** Dependence on data
# BROZ: this fails
cat("[predictCox] Dependence on data \n")
data(Melanoma, package = "riskRegression")

test_that("[predictCox] Dependence on data", {   
    Melanoma$entry <- -abs(rnorm(NROW(Melanoma), mean = 1, sd = 1))
    Melanoma2 <- Melanoma
  
    fit1 <- coxph(Surv(time, status>0)~strata(epicel)+age+strata(invasion)+sex+logthick, 
                  data = Melanoma2, x = TRUE, y = TRUE)
    GS <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  
    Melanoma2 <- 7
  
    test <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
    expect_equal(GS,test)
  
    ## with delayed entry
    ## BROZ: fails
    if (FALSE){
        Melanoma2 <- Melanoma
        fit1 <- coxph(Surv(entry ,time, status>0)~strata(epicel)+age+strata(invasion)+sex+logthick, 
                      data = Melanoma2, x = TRUE, y = TRUE)
        GS <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
  
        Melanoma2 <- 7
        
        test <- predictCox(fit1,newdata=Melanoma[1:10,],times=1000)
        expect_equal(GS,test)
    }
})


######################################################################
### test-failures.R ends here
