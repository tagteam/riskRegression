### test-ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 15 2018 (11:42) 
## Version: 
## Last-Updated: aug 31 2018 (17:12) 
##           By: Brice Ozenne
##     Update #: 34
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
                  "IPWnaive2" = c(0.15215, 0.69063, 0.53848), 
                  "IPWefficient2" = c(0.15213, 0.68416, 0.53203), 
                  "AIPWnaive" = c(0.1463, 0.75177, 0.60547), 
                  "AIPWefficient" = c(0.14629, 0.7453, 0.59901), 
                  "AIPWnaive2" = c(0.14747, 0.72998, 0.58251), 
                  "AIPWefficient2" = c(0.14746, 0.72351, 0.57605)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## ## butils::object2script(test, digit = 5)

    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.01086, 0.0194, 0.00877), 
                  "Gformula2" = c(0.04724, 0.13442, 0.12753), 
                  "IPWnaive" = c(0.05782, 0.24982, 0.26468), 
                  "IPWefficient" = c(0.0578, 0.24897, 0.2638), 
                  "IPWnaive2" = c(0.05681, 0.15178, 0.16212), 
                  "IPWefficient2" = c(0.05681, 0.1529, 0.16329), 
                  "AIPWnaive" = c(0.05662, 0.15129, 0.15937), 
                  "AIPWefficient" = c(0.05661, 0.15154, 0.15965), 
                  "AIPWnaive2" = c(0.05567, 0.14455, 0.15325), 
                  "AIPWefficient2" = c(0.05568, 0.14561, 0.15437)
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
                  "IPWnaive2" = c(0.25968, 0.48983, 0.23015), 
                  "IPWefficient2" = c(0.25948, 0.48465, 0.22517), 
                  "AIPWnaive" = c(0.25868, 0.50849, 0.24981), 
                  "AIPWefficient" = c(0.25848, 0.50331, 0.24483), 
                  "AIPWnaive2" = c(0.25982, 0.48852, 0.2287), 
                  "AIPWefficient2" = c(0.25962, 0.48334, 0.22372)
                  )
    
    expect_equal(test, M.GS, tol = 1e-4)
 
   test <- e.ateRobust$ate.se
   rownames(test) <- NULL
    
   M.GS <- cbind("Gformula" = c(0.01006, 0.01065, 6e-04), 
                 "Gformula2" = c(0.04577, 0.15642, 0.15524), 
                 "IPWnaive" = c(0.04731, 0.22229, 0.233), 
                 "IPWefficient" = c(0.04727, 0.22165, 0.23231), 
                 "IPWnaive2" = c(0.0465, 0.16437, 0.17088), 
                 "IPWefficient2" = c(0.04647, 0.16499, 0.17157), 
                 "AIPWnaive" = c(0.04636, 0.16792, 0.17302), 
                 "AIPWefficient" = c(0.04632, 0.16799, 0.1731), 
                 "AIPWnaive2" = c(0.04641, 0.16192, 0.16808), 
                 "AIPWefficient2" = c(0.0464, 0.1626, 0.16885)
                 )
   expect_equal(test, M.GS, tol = 1e-4)
   ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  "Gformula2" = c(0.28382, 0.3082, 0.02438), 
                  "IPWnaive" = c(0.25857, 0.5099, 0.25133), 
                  "IPWefficient" = c(0.25837, 0.50473, 0.24636), 
                  "IPWnaive2" = c(0.25971, 0.48983, 0.23012), 
                  "IPWefficient2" = c(0.25951, 0.48465, 0.22515), 
                  "AIPWnaive" = c(0.2587, 0.50859, 0.24988), 
                  "AIPWefficient" = c(0.2585, 0.50341, 0.24491), 
                  "AIPWnaive2" = c(0.25984, 0.48852, 0.22867), 
                  "AIPWefficient2" = c(0.25964, 0.48334, 0.2237)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.01001, 0.01058, 0.00058), 
                  "Gformula2" = c(0.04576, 0.15641, 0.15524), 
                  "IPWnaive" = c(0.04732, 0.22233, 0.23304), 
                  "IPWefficient" = c(0.04727, 0.22169, 0.23236), 
                  "IPWnaive2" = c(0.0465, 0.16443, 0.17094), 
                  "IPWefficient2" = c(0.04648, 0.16505, 0.17163), 
                  "AIPWnaive" = c(0.04636, 0.16804, 0.17314), 
                  "AIPWefficient" = c(0.04632, 0.1681, 0.17322), 
                  "AIPWnaive2" = c(0.04642, 0.16199, 0.16816), 
                  "AIPWefficient2" = c(0.0464, 0.16267, 0.16893)
                  )
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
   
})

######################################################################
### test-ateRobust.R ends here
