library(timereg)
library(mets)
library(riskRegression)

n <- 10;
d <- mets:::simCox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))

d$entry <- 0
(m1 <- phreg(Surv(entry, time,status)~X1+X2,data=d))
(m2 <- coxph(Surv(time,status)~X1+X2,data=d, x = TRUE, y = TRUE))
(m3 <- cox.aalen(Surv(time,status)~prop(X1)+prop(X2),data=d))

#### influence function ####
IC1 <- iidCox(m1)
IC2 <- iidCox(m2)

crossprod(IC1$ICbeta)
crossprod(IC2$ICbeta)
m3$robvar.gamma

#### prediction ####
predictCox(m1)
predictCox(m2)

pred1 <- predictCox(m1, newdata = d, times = 10, se = TRUE)
pred2 <- predictCox(m2, newdata = d, times = 10, se = TRUE)
pred3 <- predict(m3, newdata = d, times = 10)

data.frame(S1 = pred1$survival, S2 = pred2$survival, S3 = pred3$S0)
data.frame(S1.se = pred1$survival.se, S2.se = pred2$survival.se, S3.se = pred3$se.S0)
           

#### multiple strata
dStrata <- sampleData(1e2, outcome = "survival")
dStrata$entry <- 0
dStrata$id <- 1:NROW(dStrata)

ph.cox <- phreg(Surv(time = entry, time2 = eventtime, event = event)~strata(X1)+X2+cluster(id),data=dStrata)
# phreg(Surv(time = entry, time2 = eventtime, event = event)~strata(X1)+strata(X2)+X6+cluster(id),data=dStrata)

#### CSC ####
data(Melanoma)
fit1 <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),
            data=Melanoma, iid = TRUE)

fit2 <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),
            data=Melanoma, iid = TRUE, fitter = "phreg")
coef(fit1$models$`Cause 1`)
coef(fit2$models$`Cause 1`)

coef(fit1$models$`Cause 2`)
coef(fit2$models$`Cause 2`)

p1 <- predict(fit1, newdata = Melanoma, times = 3000, cause = 1)
p2 <- predict(fit2, newdata = Melanoma, times = 3000, cause = 1)

range(p1-p2)


#### issue with delayed entry ####
library(mets)
library(riskRegression)

n <- 100
d <- sampleData(n, outcome = "competing.risks")
d$entry <- 0

head(d)

m2 <- coxph(Surv(entry,time,event==1)~X1+X2,data=d, x = TRUE)
CSC(Hist(time = time, event = event)~ 1, data = d)
CSC(Hist(entry,time,event)~X1+X2, data = d) # error



