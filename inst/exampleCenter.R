library(survival)
library(butils.base)
package.source("riskRegression", Ccode = TRUE, RorderDescription = FALSE)


#### Baseline hazard #### 
set.seed(10)
d <- SimSurv(1e2)
nd <- SimSurv(10)
d$time <- round(d$time,1)

d$X1scaled <- scale(d$X1)
d$X2scaled <- scale(d$X2)

## Cox
CoxS <- coxph(Surv(time,status)~X1 + X2, data=d, ties="breslow")

# by default Cox center the linear predictor 
lp <- as.matrix(d[,c("X1","X2")]) %*% coef(CoxS)
range(lp - CoxS$linear.predictors)
lpS <- scale(as.matrix(d[,c("X1","X2")]), center = TRUE, scale = FALSE) %*% coef(CoxS)
range(lpS - CoxS$linear.predictors)

## Lambda
Lambda0 <- baseHaz_cpp(alltimes = CoxS$y[,"time"],
                       status = CoxS$y[,"status"],
                       eXb = exp(CoxS$linear.predictors),
                       strata = rep(0,CoxS$n),
                       se = FALSE,
                       data = matrix(0),
                       nVar = length(coef(CoxS)), 
                       nPatients = CoxS$n,
                       nStrata = 1,
                       emaxtimes = max(CoxS$y[,"time"]),
                       predtimes = numeric(0),
                       cause = 1,
                       Efron = FALSE)
range(predictCox(CoxS)$hazard - Lambda0$hazard)

#### Influence function #### 
# center.eXb: 
# center.LPdata: center the design matrix - not needed
# center.lambda0: correct computation of the baseline hazard
# center.result: corrective factor 

resRaw <- iidCox(CoxS, center.result = FALSE)
resScaled <- iidCox(CoxS, center.result = TRUE)

range(resRaw$ICbeta-resScaled$ICbeta)
range(resRaw$ICLambda0/resScaled$ICLambda0)
