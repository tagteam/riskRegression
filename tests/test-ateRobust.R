### test-ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 15 2018 (11:42) 
## Version: 
## Last-Updated: Sep 11 2018 (13:19) 
##           By: Thomas Alexander Gerds
##     Update #: 51
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

if (FALSE){
#### Survival ####

## generative model
mSimSurv <- lava::lvm()
lava::distribution(mSimSurv, ~a) <- lava::binomial.lvm(p = c(0.5))
lava::distribution(mSimSurv, ~w) <- lava::binomial.lvm(p = c(0.5))
lava::distribution(mSimSurv, "eventtime") <- lava::coxExponential.lvm(scale = 1)
lava::distribution(mSimSurv, "censtime") <- lava::coxExponential.lvm(scale = 1)
mSimSurv <- lava::eventTime(mSimSurv, time ~ min(eventtime = 1, censtime = 0), "event")
lava::regression(mSimSurv) <- eventtime ~ alpha*a + beta*w

## settings
alpha <- 0.3 
beta <- 0.4
n <- 1e3

## check bias 
set.seed(10)
dt <- as.data.table(lava::sim(mSimSurv, n = n, p = c(alpha = alpha, beta = beta)))
setkeyv(dt, c("a","w"))

    ## True value
    psi.TRUE <- c("risk.0" = (1-exp(-1*exp(alpha*0+beta*0)))*0.5 + (1-exp(-1*exp(alpha*0+beta*1)))*0.5,
                  "risk.1" = (1-exp(-1*exp(alpha*1+beta*0)))*0.5 + (1-exp(-1*exp(alpha*1+beta*1)))*0.5)
    psi.TRUE

    ## Approximate true value
    dt[,.(risk = mean(eventtime<1)),by = c("a")]

    ## Estimated using stratified Cox model
    res1 <- ateRobust(data = dt, type = "survival",
                      formula.event = Surv(time, event) ~ strata(a,w),
                      formula.censor = Surv(time, event==0) ~ strata(a,w),
                      formula.treatment = a ~ w,
                      times = 1,
                      product.limit = FALSE)
    print(res1, augment.cens = TRUE, nuisance.iid = TRUE)
    print(res1, augment.cens = TRUE, nuisance.iid = FALSE)
    print(res1, augment.cens = FALSE, nuisance.iid = TRUE)
    print(res1, augment.cens = FALSE, nuisance.iid = FALSE)
}
if(FALSE){
    e.cox <- coxph(Surv(time, event) ~ strata(a,w), data = dt, x = TRUE, y = TRUE)
    ate(e.cox, data = dt, treatment = "a", times = 1)
    rbind(estimate = res1$ate.value[,"Gformula"],
          se = res1$ate.se[,"Gformula"])

## Estimate using Cox model
res2 <- ateRobust(data = dt, type = "survival",
            formula.event = Surv(time, event) ~ a + w,
            formula.censor = Surv(time, event==0) ~ a + w,
            formula.treatment = a ~ w,
            times = 1,
            product.limit = FALSE,
            nuisance.iid = TRUE)

print(res2, augment.cens = TRUE, nuisance.iid = TRUE)
print(res2, augment.cens = TRUE, nuisance.iid = FALSE)
print(res2, augment.cens = FALSE, nuisance.iid = TRUE)
print(res2, augment.cens = FALSE, nuisance.iid = FALSE)

#### Competing risks ####
set.seed(10)
n <- 1e3

## simulate data
alphaE.X <- 2
alphaCR.X <- 1
alphaE.Y <- 3
alphaCR.Y <- 2
set.seed(10)
df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X),
                       time2 = rexp(n, rate = alphaCR.X), group = "1"),
            data.frame(time1 = rexp(n, rate = alphaE.Y),
                       time2 = rexp(n, rate = alphaCR.Y), group = "2"))
df$time <- pmin(df$time1,df$time2) ## first event
df$event <- (df$time2<df$time1)+1 ## type of event
df$eventC <- df$event
df$eventC[rbinom(n, size = 1, prob = 0.2)==1] <- 0

## true value
tau <- 1
c(CIF.X = alphaE.X/(alphaE.X+alphaCR.X)*(1-exp(-(alphaE.X+alphaCR.X)*(tau))),
  CIF.Y = alphaE.Y/(alphaE.Y+alphaCR.Y)*(1-exp(-(alphaE.Y+alphaCR.Y)*(tau))))
## estimating using a CSC (no censoring)
tapply(df$time,df$group,max)

    resCR <- ateRobust(data = df, type = "competing.risks",
                       formula.event = Hist(time, event) ~ group, ## strata(group),
                       formula.censor = Surv(time, event==0) ~ group,## strata(group),
                       formula.treatment = group ~ 1,
                       times = tau,
                       nuisance.iid = FALSE,
                       product.limit = FALSE,
                       cause = 1)
    resCR
    ## estimating using a CSC (censoring)
    resCRC <- ateRobust(data = df, type = "competing.risks",
                        formula.event = Hist(time, eventC) ~  strata(group), ## group, ##
                        formula.censor = Surv(time, eventC==0) ~  strata(group), ## group,##
                        formula.treatment = group ~ 1,
                        times = tau,
                        nuisance.iid = FALSE,
                        product.limit = FALSE,
                        cause = 1)
    print(resCRC, augment.cens = TRUE)
    print(resCRC, augment.cens = FALSE)
}

## * survival case

    if (FALSE){
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
                             type = "survival")
    
    ## e.ateRobust$ate.value
    expect_equal(as.double(e.ateRobust$ate.value[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk],-e.ate$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobust$ate.se[,"Gformula"]),
                 c(e.ate$meanRisk[,meanRisk.se],e.ate$riskComparison[,diff.se]))

    a <- c("a","b","c","d","e","f","g","h","i","j","k"
          ,"l","m","n","o","p","q","r","s","t","u","v",
           "w","x","y","z")


    ## ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.19447, 0.54518, 0.35071), 
                  "IPTW.IPCW" = c(0.15098, 0.71242, 0.56145), 
                  "IPTW.AIPCW" = c(0.15097, 0.70596, 0.55499), 
                  "AIPTW.IPCW" = c(0.1463, 0.75177, 0.60547), 
                  "AIPTW.AIPCW" = c(0.14629, 0.7453, 0.59901))
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)

    test <- e.ateRobust$ate.se
    rownames(test) <- NULL
    M.GS <- cbind("Gformula" = c(0.04724, 0.13442, 0.12753), 
                  "IPTW.IPCW" = c(0.05782, 0.24982, 0.26468), 
                  "IPTW.AIPCW" = c(0.0578, 0.24897, 0.2638), 
                  "IPTW.IPCW" = c(0.05683, 0.14767, 0.15823), 
                  "AIPTW.IPCW" = c(0.05662, 0.15129, 0.15937), 
                  "AIPTW.AIPCW" = c(0.05661, 0.15154, 0.15965))
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
                             product.limit = FALSE)
    e.ateRobustPL <- ateRobust(data = dtS, times = 3, cause = 1,
                             formula.event = Hist(time,event) ~ X1 + X2 + X3,
                             formula.censor = Surv(time,event==0) ~ X1,
                             formula.treatment = X1 ~ 1,
                             type = "competing.risks",
                             product.limit = TRUE)

    expect_equal(as.double(e.ateRobustPL$ate.value[,"Gformula"]),
                 c(e.atePL$meanRisk[,meanRisk],-e.atePL$riskComparison[,diff]))
    expect_equal(as.double(e.ateRobustPL$ate.se[,"Gformula"]),
                 c(e.atePL$meanRisk[,meanRisk.se],e.atePL$riskComparison[,diff.se]))

    ## check values
    test <- e.ateRobust$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28427, 0.30878, 0.02452), 
                  "IPTW.IPCW" = c(0.25855, 0.5098, 0.25125), 
                  "IPTW.AIPCW" = c(0.25835, 0.50462, 0.24627), 
                  "AIPTW.IPCW" = c(0.25868, 0.50849, 0.24981), 
                  "AIPTW.AIPCW" = c(0.25848, 0.50331, 0.24483))
    
    expect_equal(test, M.GS, tol = 1e-4)
   ## butils::object2script(test, digit = 5) 
   test <- e.ateRobust$ate.se
   rownames(test) <- NULL
    
   M.GS <- cbind("Gformula" = c(0.01006, 0.01065, 6e-04), 
                 "IPTW.IPCW" = c(0.04731, 0.22229, 0.233), 
                 "IPTW.AIPCW" = c(0.04727, 0.22165, 0.23231), 
                 "AIPTW.IPCW" = c(0.04636, 0.16792, 0.17302), 
                 "AIPTW.AIPCW" = c(0.04632, 0.16799, 0.1731))
   expect_equal(test, M.GS, tol = 1e-4)
   ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.value
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.28382, 0.3082, 0.02438), 
                  "IPTW.IPCW" = c(0.25857, 0.5099, 0.25133), 
                  "IPTW.AIPCW" = c(0.25837, 0.50473, 0.24636), 
                  "AIPTW.IPCW" = c(0.2587, 0.50859, 0.24988), 
                  "AIPTW.AIPCW" = c(0.2585, 0.50341, 0.24491))
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
    
    test <- e.ateRobustPL$ate.se
    rownames(test) <- NULL
    
    M.GS <- cbind("Gformula" = c(0.01001, 0.01058, 0.00058), 
                  "IPTW.IPCW" = c(0.04732, 0.22233, 0.23304), 
                  "IPTW.AIPCW" = c(0.04727, 0.22169, 0.23236), 
                  "AIPTW.IPCW" = c(0.04636, 0.16804, 0.17314), 
                  "AIPTW.AIPCW" = c(0.04632, 0.1681, 0.17322))
    expect_equal(test, M.GS, tol = 1e-4)
    ## butils::object2script(test, digit = 5)
   
})
}
######################################################################
### test-ateRobust.R ends here
