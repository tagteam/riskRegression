### simMelanoma.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 27 2017 (09:52) 
## Version: 
## last-updated: Feb 27 2017 (12:54) 
##           By: Thomas Alexander Gerds
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Simulate data alike the Melanoma data
##'
##' This is based on the functionality of \code{library(lava)}.
##' @title Simulate data alike the Melanoma data
##' @param n sample size
##' @return data table of size n
##' @examples
##' set.seed(71)
##' simMelanoma(3)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
simMelanoma <- function(n){
    library(lava)
    library(survival)
    library(riskRegression)
    data(Melanoma,package="riskRegression")
    Melanoma$ici1 <- 1*(Melanoma$ici=="1")
    Melanoma$ici2 <- 1*(Melanoma$ici=="2")
    Melanoma$ici3 <- 1*(Melanoma$ici=="3")
    pici <- prop.table(table(Melanoma$ici))
    pulcer <- prop.table(table(Melanoma$ulcer))
    pepicel <- prop.table(table(Melanoma$epicel))
    psex <- prop.table(table(Melanoma$sex))
    l0 <- glm(sex~ulcer+epicel+ici1+ici2+ici3,data=Melanoma,family="binomial")
    l1 <- glm(logthick~age+sex+ulcer+epicel+ici1+ici2+ici3,data=Melanoma)
    l2 <- glm(age~sex,data=Melanoma)
    G <- survreg(Surv(time,status==0)~1,data=Melanoma)
    F1 <- survreg(Surv(time,status==1)~logthick+age+sex+ulcer+epicel+ici1+ici2+ici3,data=Melanoma)
    F2 <- survreg(Surv(time,status==2)~age+sex,data=Melanoma)
    m <- lvm(~ulcer+epicel+sex+age+logthick)
    ## m <- categorical(m,~ici, K=4,p=pici[-3])
    distribution(m,~ici1) <- binomial.lvm(p=pici[2])
    distribution(m,~ici2) <- binomial.lvm(p=pici[3])
    distribution(m,~ici3) <- binomial.lvm(p=pici[4])
    distribution(m,~sex) <- binomial.lvm(p=psex[2])
    distribution(m,~ulcer) <- binomial.lvm(p=pulcer[2])
    distribution(m,~epicel) <- binomial.lvm(p=pepicel[2])
    distribution(m,~t0) <- coxWeibull.lvm(scale=exp(-2-G$coef["(Intercept)"]/G$scale),shape=1/G$scale)
    distribution(m,~t1) <- coxWeibull.lvm(scale=exp(-F1$coef["(Intercept)"]/F1$scale),shape=1/F1$scale)
    distribution(m,~t2) <- coxWeibull.lvm(scale=exp(3+-F2$coef["(Intercept)"]/F2$scale),shape=1/F2$scale)
    transform(m,thick~logthick) <- function(x){exp(x)}
    regression(m,sex~ulcer+epicel+ici1+ici2+ici3) <- coef(l0)[-1]
    regression(m,logthick~age+sex+ulcer+epicel+ici1+ici2+ici3) <- coef(l1)[-1]
    regression(m,age~sex) <- coef(l2)[-1]
    regression(m,t1~logthick+age+sex+ulcer+epicel+ici1+ici2+ici3) <- -coef(F1)[-1]/F1$scale
    regression(m,t2~age+sex) <- -coef(F2)[-1]/F2$scale
    m <- eventTime(m,time~min(t0=0,t1=1,t2=2),"status")
    simdat <- sim(m,n)
    setDT(simdat)
    simdat[,ici:=factor("0",levels=0:3)]
    simdat[ici1==1,ici:="1"]
    simdat[ici2==1,ici:="2"]
    simdat[ici3==1,ici:="3"]
    simdat[,c("ici1","ici2","ici3","t0","t1","t2"):=NULL]
    simdat$ici <- factor(simdat$ici)
    simdat$ulcer <- factor(simdat$ulcer,levels=c("0","1"),labels=c("not present","present"))
    simdat$epicel <- factor(simdat$epicel,levels=c("0","1"),labels=c("not present","present"))
    simdat$sex <- factor(simdat$sex,levels=c("0","1"),labels=c("Female","Male"))
    simdat
}


#----------------------------------------------------------------------
### simMelanoma.R ends here
