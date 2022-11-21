### simPBC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  1 2022 (18:03) 
## Version: 
## Last-Updated: Nov  1 2022 (18:20) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' simulating data alike the pbc data
##'
##' using lava to synthesize data
##' @title bla
##' @param n Sample size
##' @return bla 
##' @seealso bla
##' @export simPBC
##' @examples 
##' library(survival)
##' library(lava)
##' # simulate data alike pbc data
##' set.seed(98)
##' d=simPBC(847)
##' d$protimegrp1 <- d$protimegrp=="10-11"
##' d$protimegrp2 <- d$protimegrp==">11"
##' d$sex <- factor(d$sex,levels=0:1,labels=c("m","f"))
##' sF1 <- survreg(Surv(time,status==1)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=d)
##' coxF1 <- coxph(Surv(time,status==1)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=d)
##' # load real pbc data
##' data(pbc,package="survival")
##' pbc <- na.omit(pbc[,c("time","status","age","sex","stage","bili","protime","trt")])
##' pbc$stage <- factor(pbc$stage)
##' levels(pbc$stage) <- list("1/2"=c(1,2),"3"=3,"4"=4)
##' pbc$logbili <- log(pbc$bili)
##' pbc$logprotime <- log(pbc$protime)
##' pbc$stage3 <- 1*(pbc$stage=="3")
##' pbc$stage4 <- 1*(pbc$stage=="4")
##' pbc$protimegrp <- cut(pbc$protime,c(-Inf,10,11,Inf),labels=c("<=10","10-11",">11"))
##' pbc$protimegrp1 <- pbc$protimegrp=="10-11"
##' pbc$protimegrp2 <- pbc$protimegrp==">11"
##' F1 <- survival::survreg(Surv(time,status==1)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
##' F2 <- survival::survreg(Surv(time,status==2)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
##' sF2 <- survreg(Surv(time,status==2)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=d)
##' G <- survreg(Surv(time,status==0)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
##' sG <- survreg(Surv(time,status==0)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=d)
##' # compare fits in real and simulated pbc data
##' cbind(coef(F1),coef(sF1))
##' cbind(coef(F2),coef(sF2))
##' cbind(coef(G),coef(sG))
##' cbind(coef(glm(protimegrp1~age+sex+logbili,data=pbc,family="binomial")),coef(glm(protimegrp1~age+sex+logbili,data=d,family="binomial")))
##' cbind(coef(lm(logbili~age+sex,data=pbc)),coef(lm(logbili~age+sex,data=d)))
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
simPBC <- function(n){
    requireNamespace("survival")
    requireNamespace("lava")
    ## library(survival)
    ## library(prodlim)
    ## library(lava)
    data(pbc,package = "survival")
    pbc$stage <- factor(pbc$stage)
    levels(pbc$stage) <- list("1/2"=c(1,2),"3"=3,"4"=4)
    pbc <- na.omit(pbc[,c("time","status","age","sex","stage","bili","protime","trt")])
    pbc$logbili <- log(pbc$bili)
    pbc$logprotime <- log(pbc$protime)
    pbc$stage3 <- 1*(pbc$stage=="3")
    pbc$stage4 <- 1*(pbc$stage=="4")
    pbc$protimegrp <- cut(pbc$protime,c(-Inf,10,11,Inf),labels=c("<=10","10-11",">11"))
    pbc$protimegrp1 <- pbc$protimegrp=="10-11"
    pbc$protimegrp2 <- pbc$protimegrp==">11"
    ## f1 <- coxph(Surv(time,status==1)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4  ,data=pbc)
    ## f2 <- coxph(Surv(time,status==2)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4  ,data=pbc)
    ## g <- coxph(Surv(time,status==0)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4   ,data=pbc)
    ## parametric survival
    F1 <- survival::survreg(Surv(time,status==1)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
    ## F1a <- survival::survreg(Surv(time,status==1)~sex+age+logbili+protimegrp+stage,data=pbc)
    F2 <- survival::survreg(Surv(time,status==2)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
    G <- survreg(Surv(time,status==0)~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4,data=pbc)
    ## G <- survival::survreg(Surv(time,status==0)~1,data=pbc)
    ## logistic lava::regression
    M <- lava::lvm()
    lava::distribution(M,~sex) <- lava::binomial.lvm(p=mean(pbc$sex=="f"))
    lava::distribution(M,~age) <- lava::normal.lvm(mean=50,sd=10)
    lava::distribution(M,~trt) <- lava::binomial.lvm(p=.5)
    lava::distribution(M,~logbili) <- lava::normal.lvm(mean=mean(pbc$logbili),sd=sd(pbc$logbili))
    lava::distribution(M,~logprotime) <- lava::normal.lvm(mean=mean(pbc$logprotime),sd=sd(pbc$logprotime))
    lava::distribution(M,~time.transplant) <- lava::coxWeibull.lvm(scale=exp(-F1$coef["(Intercept)"]/F1$scale),shape=1/F1$scale)
    lava::distribution(M,~time.death) <- lava::coxWeibull.lvm(scale=exp(-F2$coef["(Intercept)"]/F2$scale),shape=1/F2$scale)
    lava::distribution(M,~time.cens) <- lava::coxWeibull.lvm(scale=exp(-G$coef["(Intercept)"]/G$scale),shape=1/G$scale)
    ## bili coefficients
    lava::regression(M,logbili~age+sex) <- coef(lm(logbili~age+sex,data=pbc))[-1]
    ## protime coefficients
    lava::regression(M,protimegrp1~age+sex+logbili) <- coef(glm(protimegrp1~age+sex+logbili,data=pbc,family="binomial"))[-1]
    lava::regression(M,protimegrp2~age+sex+logbili) <- coef(glm(protimegrp2~age+sex+logbili,data=pbc,family="binomial"))[-1]
    ## stage coefficients
    lava::regression(M,stage3~age+sex+protimegrp1+protimegrp2+logbili) <- coef(glm(stage3~age+sex+protimegrp1+protimegrp2+logbili,data=pbc,family="binomial"))[-1]
    lava::regression(M,stage4~age+sex+protimegrp1+protimegrp2+logbili) <- coef(glm(stage4~age+sex+protimegrp1+protimegrp2+logbili,data=pbc,family="binomial"))[-1]
    ## survival coefficients
    lava::regression(M,time.transplant~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4) <- -coef(F1)[-1]/F1$scale
    lava::regression(M,time.death~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4) <- -coef(F2)[-1]/F2$scale
    lava::regression(M,time.cens~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4) <- -coef(G)[-1]/G$scale
    M <- lava::categorical(M,~stage,labels=c("1/2","3","4"),K=3,p=c(mean(pbc$stage=="3"),mean(pbc$stage=="4")))
    lava::transform(M,stage3~stage) <- function(x){1*(x==1)}
    lava::transform(M,stage4~stage) <- function(x){1*(x==2)}
    lava::transform(M,protime~logprotime) <- function(x){exp(x)}
    lava::transform(M,protimegrp~protime) <- function(x){cut(x[,"protime"], c(-Inf,10,11,Inf), labels=c("<=10","10-11",">11"))}
    lava::transform(M,protimegrp1~protimegrp) <- function(x){1*(x=="10-11")}
    lava::transform(M,protimegrp2~protimegrp) <- function(x){1*(x==">11")}
    M <- lava::eventTime(M,time~min(time.cens=0,time.transplant=1,time.death=2),"status")
    d <- lava::sim(M,n)
    d
}


######################################################################
### simPBC.R ends here
