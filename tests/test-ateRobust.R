### test-ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 15 2018 (11:42) 
## Version: 
## Last-Updated: aug 20 2018 (22:01) 
##           By: Brice Ozenne
##     Update #: 24
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

    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk],-e.ate$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula2"]),
                 as.double(e.ateRobust$ate.value[,"Gformula"]))
    expect_equal(as.double(e.ateRobust$ate.se[,"Gformula2"]),
                 c(e.ate$meanRisk[,meanRisk.se],e.ate$riskComparison[,diff.se]))


    ## ## check values
    ## test <- e.ateRobust$ate.value
    ## rownames(test) <- NULL
    ## M.GS <- cbind("Gformula" = c(0.19447, 0.54518, 0.35071), 
    ##               "Gformula2" = c(0.19447, 0.54518, 0.35071), 
    ##               "IPWnaive" = c(0.15098, 0.71242, 0.56145), 
    ##               "IPWefficient" = c(0.15061, 0.70431, 0.5537), 
    ##               "IPWnaive2" = c(0.15215, 0.69063, 0.53848), 
    ##               "IPWefficient2" = c(0.15178, 0.68252, 0.53074), 
    ##               "AIPWnaive" = c(0.14726, 0.73904, 0.59178), 
    ##               "AIPWefficient" = c(0.14689, 0.73092, 0.58403), 
    ##               "AIPWnaive2" = c(0.14686, 0.73285, 0.58599), 
    ##               "AIPWefficient2" = c(0.14649, 0.72473, 0.57824)
    ##               )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## ## butils::object2script(test, digit = 5)

    ## test <- e.ateRobust$ate.se
    ## rownames(test) <- NULL
    ## M.GS <- cbind("Gformula" = c(0.01086, 0.0194, 0.00877), 
    ##               "Gformula2" = c(0.04724, 0.13442, 0.12753), 
    ##               "IPWnaive" = c(0.05782, 0.24982, 0.26468), 
    ##               "IPWefficient" = c(0.05781, 0.24876, 0.26357), 
    ##               "IPWnaive2" = c(0.05681, 0.15178, 0.16212), 
    ##               "IPWefficient2" = c(0.05682, 0.15319, 0.16357), 
    ##               "AIPWnaive" = c(0.05661, 0.1523, 0.16035), 
    ##               "AIPWefficient" = c(0.05663, 0.15274, 0.16089), 
    ##               "AIPWnaive2" = c(0.05569, 0.14401, 0.15266), 
    ##               "AIPWefficient2" = c(0.0557, 0.14534, 0.15403)
    ##               )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## ## butils::object2script(test, digit = 5)
})


## * competing risk case
set.seed(10)
n <- 1e2
dtS <- sampleData(n,outcome="competing.risks")
## dtS[,min(time),by = event]

## ** agreement with ate
e.CSC <- CSC(Hist(time, event) ~ X1 + X2 + X3,
             data = dtS, surv.type = "hazard")

test_that("Agreement ate-ateRobust (survival)",{
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
    
    ## M.GS <- cbind("Gformula" = c(0.28427, 0.30878, 0.02452), 
                  ## "Gformula2" = c(0.28427, 0.30878, 0.02452), 
                  ## "IPWnaive" = c(0.25859, 0.50988, 0.25129), 
                  ## "IPWefficient" = c(0.25823, 0.50809, 0.24986), 
                  ## "IPWnaive2" = c(0.23861, 0.50875, 0.27013), 
                  ## "IPWefficient2" = c(0.23826, 0.50696, 0.2687), 
                  ## "AIPWnaive" = c(0.2598, 0.49796, 0.23816), 
                  ## "AIPWefficient" = c(0.25944, 0.49618, 0.23673), 
                  ## "AIPWnaive2" = c(0.25059, 0.49826, 0.24768), 
                  ## "AIPWefficient2" = c(0.25023, 0.49648, 0.24625)
                  ## )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    
    ## M.GS <- cbind("Gformula" = c(0.01006, 0.01065, 6e-04), 
                  ## "Gformula2" = c(0.04577, 0.15642, 0.15524), 
                  ## "IPWnaive" = c(0.04732, 0.22232, 0.23303), 
                  ## "IPWefficient" = c(0.04731, 0.22206, 0.23275), 
                  ## "IPWnaive2" = c(0.17079, 0.21647, 0.34982), 
                  ## "IPWefficient2" = c(0.17075, 0.21622, 0.34939), 
                  ## "AIPWnaive" = c(0.04636, 0.16832, 0.17344), 
                  ## "AIPWefficient" = c(0.04636, 0.16833, 0.17347), 
                  ## "AIPWnaive2" = c(0.08055, 0.16874, 0.20275), 
                  ## "AIPWefficient2" = c(0.08052, 0.16874, 0.20259)
                  ## )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    ## M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  ## "Gformula2" = c(0.28382, 0.3082, 0.02438), 
                  ## "IPWnaive" = c(0.25855, 0.50987, 0.25132), 
                  ## "IPWefficient" = c(0.2582, 0.50809, 0.24989), 
                  ## "IPWnaive2" = c(0.23849, 0.50873, 0.27025), 
                  ## "IPWefficient2" = c(0.23813, 0.50695, 0.26882), 
                  ## "AIPWnaive" = c(0.25977, 0.49791, 0.23813), 
                  ## "AIPWefficient" = c(0.25942, 0.49612, 0.23671), 
                  ## "AIPWnaive2" = c(0.2505, 0.49821, 0.24771), 
                  ## "AIPWefficient2" = c(0.25014, 0.49642, 0.24628)
                  ## )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    ## M.GS <- cbind("Gformula" = c(0.01001, 0.01058, 0.00058), 
                  ## "Gformula2" = c(0.04576, 0.15641, 0.15524), 
                  ## "IPWnaive" = c(0.04731, 0.22232, 0.23303), 
                  ## "IPWefficient" = c(0.0473, 0.22206, 0.23275), 
                  ## "IPWnaive2" = c(0.17081, 0.21647, 0.34982), 
                  ## "IPWefficient2" = c(0.17078, 0.21622, 0.34939), 
                  ## "AIPWnaive" = c(0.04636, 0.1684, 0.17352), 
                  ## "AIPWefficient" = c(0.04635, 0.16841, 0.17355), 
                  ## "AIPWnaive2" = c(0.08068, 0.16881, 0.20295), 
                  ## "AIPWefficient2" = c(0.08065, 0.16881, 0.20278)
                  ## )
    ## expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    

})

######################################################################
### test-ateRobust.R ends here
