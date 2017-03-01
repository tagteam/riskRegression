### simMelanoma.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 27 2017 (09:52) 
## Version: 
## last-updated: Feb 28 2017 (20:09) 
##           By: Thomas Alexander Gerds
##     Update #: 26
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
    ici=ici1=ici2=ici3=NULL
    ## data(Melanoma,package="riskRegression")
    ## Melanoma$ici1 <- 1*(Melanoma$ici=="1")
    ## Melanoma$ici2 <- 1*(Melanoma$ici=="2")
    ## Melanoma$ici3 <- 1*(Melanoma$ici=="3")
    ## pici <- prop.table(table(Melanoma$ici))
    ## pulcer <- prop.table(table(Melanoma$ulcer))
    ## pepicel <- prop.table(table(Melanoma$epicel))
    ## psex <- prop.table(table(Melanoma$sex))
    ## l0 <- glm(sex~ulcer+epicel+ici1+ici2+ici3,data=Melanoma,family="binomial")
    ## l1 <- glm(logthick~age+sex+ulcer+epicel+ici1+ici2+ici3,data=Melanoma)
    ## l2 <- glm(age~sex,data=Melanoma)
    ## G <- survival::survreg(Surv(time,status==0)~1,data=Melanoma)
    ## F1 <- survival::survreg(Surv(time,status==1)~logthick+age+sex+ulcer+epicel+ici1+ici2+ici3,data=Melanoma)
    ## F2 <- survival::survreg(Surv(time,status==2)~age+sex,data=Melanoma)
    m <- lava::lvm(~ulcer+epicel+sex+age+logthick)
    ## m <- categorical(m,~ici, K=4,p=pici[-3])
    lava::distribution(m,~ici1) <- lava::binomial.lvm(p=0.29)
    lava::distribution(m,~ici2) <- lava::binomial.lvm(p=0.52)
    lava::distribution(m,~ici3) <- lava::binomial.lvm(p=0.11)
    lava::distribution(m,~sex) <- lava::binomial.lvm(p=0.39)
    lava::distribution(m,~ulcer) <- lava::binomial.lvm(p=0.44)
    lava::distribution(m,~epicel) <- lava::binomial.lvm(p=0.44)
    lava::distribution(m,~t0) <- lava::coxWeibull.lvm(scale=exp(-2-8/0.33),shape=1/0.33)
    lava::distribution(m,~t1) <- lava::coxWeibull.lvm(scale=exp(-11.5/0.8),shape=1/0.8)
    lava::distribution(m,~t2) <- lava::coxWeibull.lvm(scale=exp(3+-16/1.3),shape=1/1.3)
    transform(m,thick~logthick) <- function(x){exp(x)}
    lava::regression(m,sex~ulcer+epicel+ici1+ici2+ici3) <- c(0.83,0.67,-0.28,0.03,-0.28)
    lava::regression(m,logthick~age+sex+ulcer+epicel+ici1+ici2+ici3) <- c(0.01,0.23,0.9,-0.35,0.56,0.73,0.83)
    lava::regression(m,age~sex) <- 2.34
    lava::regression(m,t1~logthick+age+sex+ulcer+epicel+ici1+ici2+ici3) <- -c(-0.29,-0.02,-0.34,-0.77,0.53,-1.39,-1.32,-1.7)/0.8
    lava::regression(m,t2~age+sex) <- -c(-0.08,-0.39)/1.3
    m <- lava::eventTime(m,time~min(t0=0,t1=1,t2=2),"status")
    simdat <- lava::sim(m,n)
    data.table::setDT(simdat)
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
