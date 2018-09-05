### test-ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 15 2018 (11:42) 
## Version: 
## Last-Updated: sep  3 2018 (09:35) 
##           By: Brice Ozenne
##     Update #: 39
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(riskRegression)
library(survival)
library(testthat)
context("Ate robust checks")

## * survival case

## ** agreement with ate

set.seed(10)
n <- 5e1
dtS <- sampleData(n,outcome="survival")
e.cox <- coxph(Surv(time, event) ~ X1 + X2 + X3,
               data = dtS,
               x = TRUE)

test_that("Agreement ate-ateRobust (survival)",{
    e.ate <- ate(e.cox, treatment = "X1", times = 3, data = dtS, se = TRUE)
    e.ateRobust <- ateRobust(data = dtS, times = 3,
                             formula.event = Surv(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "survival",
                             nuisance.iid = TRUE)
    
    ## e.ateRobust$ate.value
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk],-e.ate$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula2"]),
                 as.double(e.ateRobust$ate.value[,"Gformula"]))
    expect_equal(as.double(e.ateRobust$ate.se[,"Gformula2"]),
                 c(e.ate$meanRisk[,meanRisk.se],e.ate$riskComparison[,diff.se]))


    ## ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.19447, 0.54518, 0.35071), 
                  "Gformula2" = c(0.19447, 0.54518, 0.35071), 
                  "IPWnaive" = c(0.15098, 0.71242, 0.56145), 
                  "IPWefficient" = c(0.15097, 0.70596, 0.55499), 
                  "IPWnaive2" = c(0.15098, 0.71242, 0.56145), 
                  "IPWefficient2" = c(0.15097, 0.70596, 0.55499), 
                  "AIPWnaive" = c(0.1463, 0.75177, 0.60547), 
                  "AIPWefficient" = c(0.14629, 0.7453, 0.59901), 
                  "AIPWnaive2" = c(0.1463, 0.75177, 0.60547), 
                  "AIPWefficient2" = c(0.14629, 0.7453, 0.59901)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)

    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.01086, 0.0194, 0.00877), 
                  "Gformula2" = c(0.04724, 0.13442, 0.12753), 
                  "IPWnaive" = c(0.05782, 0.24982, 0.26468), 
                  "IPWefficient" = c(0.0578, 0.24897, 0.2638), 
                  "IPWnaive2" = c(0.05683, 0.14767, 0.15823), 
                  "IPWefficient2" = c(0.05681, 0.14873, 0.15921), 
                  "AIPWnaive" = c(0.05662, 0.15129, 0.15937), 
                  "AIPWefficient" = c(0.05661, 0.15154, 0.15965), 
                  "AIPWnaive2" = c(0.05571, 0.14053, 0.14949), 
                  "AIPWefficient2" = c(0.0557, 0.14153, 0.15041)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
})


## * competing risk case
set.seed(10)
n <- 1e2
dtS <- sampleData(n,outcome="competing.risks")
## dtS[,min(time),by = event]

## ** agreement with ate
e.CSC <- CSC(Hist(time, event) ~ X1 + X2 + X3,
             data = dtS, surv.type = "hazard")

test_that("Agreement ate-ateRobust (competing.risks)",{
    ## NOT POSSIBLE: predictRisk does not recognise the argument product.limit
    ## e.ate <- ate(e.CSC, treatment = "X1", times = 3, data = dtS, cause = 1,
                 ## se = TRUE, product.limit = FALSE) 
    e.atePL <- ate(e.CSC, treatment = "X1", times = 3, data = dtS, cause = 1,
                   se = TRUE)

    ## data0 <- copy(dtS)
    ## data0[, X1 := factor("0", level = levels(dtS$X1))]
    ## data1 <- copy(dtS)
    ## data1[, X1 := factor("1", level = levels(dtS$X1))]
    ## mean(predict(e.CSC, cause = 1, newdata = data0, times = 3, product.limit = TRUE)$absRisk)
    
    e.ateRobust <- ateRobust(data = dtS, times = 3, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             nuisance.iid = TRUE,
                             product.limit = FALSE)
    e.ateRobustPL <- ateRobust(data = dtS, times = 3, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             nuisance.iid = TRUE,
                             product.limit = TRUE)

    expect_equal(as.double(e.ateRobustPL$ate.value[,"Gformula"]),
                 c(e.atePL$meanRisk[,meanRisk],-e.atePL$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula2"]),
                 as.double(e.ateRobust$ate.value[,"Gformula"]))
    expect_equal(as.double(e.ateRobust$atePL.value[,"Gformula2"]),
                 as.double(e.ateRobust$atePL.value[,"Gformula"]))

    expect_equal(as.double(e.ateRobustPL$ate.se[,"Gformula2"]),
                 c(e.atePL$meanRisk[,meanRisk.se],e.atePL$riskComparison[,diff.se]))

    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28427, 0.30878, 0.02452), 
                  "Gformula2" = c(0.28427, 0.30878, 0.02452), 
                  "IPWnaive" = c(0.25855, 0.5098, 0.25125), 
                  "IPWefficient" = c(0.25835, 0.50462, 0.24627), 
                  "IPWnaive2" = c(0.25855, 0.5098, 0.25125), 
                  "IPWefficient2" = c(0.25835, 0.50462, 0.24627), 
                  "AIPWnaive" = c(0.25868, 0.50849, 0.24981), 
                  "AIPWefficient" = c(0.25848, 0.50331, 0.24483), 
                  "AIPWnaive2" = c(0.25868, 0.50849, 0.24981), 
                  "AIPWefficient2" = c(0.25848, 0.50331, 0.24483)
                  )
    
    expect_equal(test, M.GS, tol = 1e-4)
   ## butils::object2script(test, digit = 5) 
   test <- e.ateRobust$ate.se
   rownames(test) <- NULL
    
   M.GS <- cbind("Gformula" = c(0.01006, 0.01065, 6e-04), 
                 "Gformula2" = c(0.04577, 0.15642, 0.15524), 
                 "IPWnaive" = c(0.04731, 0.22229, 0.233), 
                 "IPWefficient" = c(0.04727, 0.22165, 0.23231), 
                 "IPWnaive2" = c(0.04652, 0.16131, 0.16789), 
                 "IPWefficient2" = c(0.04648, 0.1619, 0.16844), 
                 "AIPWnaive" = c(0.04636, 0.16792, 0.17302), 
                 "AIPWefficient" = c(0.04632, 0.16799, 0.1731), 
                 "AIPWnaive2" = c(0.04644, 0.15882, 0.16504), 
                 "AIPWefficient2" = c(0.0464, 0.15946, 0.16567)
                 )
   expect_equal(test, M.GS, tol = 1e-4)
   ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  "Gformula2" = c(0.28382, 0.3082, 0.02438), 
                  "IPWnaive" = c(0.25857, 0.5099, 0.25133), 
                  "IPWefficient" = c(0.25837, 0.50473, 0.24636), 
                  "IPWnaive2" = c(0.25857, 0.5099, 0.25133), 
                  "IPWefficient2" = c(0.25837, 0.50473, 0.24636), 
                  "AIPWnaive" = c(0.2587, 0.50859, 0.24988), 
                  "AIPWefficient" = c(0.2585, 0.50341, 0.24491), 
                  "AIPWnaive2" = c(0.2587, 0.50859, 0.24988), 
                  "AIPWefficient2" = c(0.2585, 0.50341, 0.24491)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.01001, 0.01058, 0.00058), 
                  "Gformula2" = c(0.04576, 0.15641, 0.15524), 
                  "IPWnaive" = c(0.04732, 0.22233, 0.23304), 
                  "IPWefficient" = c(0.04727, 0.22169, 0.23236), 
                  "IPWnaive2" = c(0.04653, 0.16134, 0.16792), 
                  "IPWefficient2" = c(0.04648, 0.16193, 0.16847), 
                  "AIPWnaive" = c(0.04636, 0.16804, 0.17314), 
                  "AIPWefficient" = c(0.04632, 0.1681, 0.17322), 
                  "AIPWnaive2" = c(0.04644, 0.15886, 0.16509), 
                  "AIPWefficient2" = c(0.0464, 0.15951, 0.16572)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
   
})

######################################################################
### test-ateRobust.R ends here
