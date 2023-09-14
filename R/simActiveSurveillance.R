### simActiveSurveillance.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 27 2017 (09:52) 
## Version: 
## last-updated: May 26 2019 (15:38) 
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
##' Simulate data of a hypothetical active surveillance prostate cancer study
##'
##' This is based on the functionality of \code{library(lava)}.
##' @title Simulate data of a hypothetical active surveillance prostate cancer study
##' @param n sample size
##' @return data table of size n
##' @examples
##' set.seed(71)
##' simActiveSurveillance(3)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
simActiveSurveillance <- function(n){
    ## setwd("~/research/Methods/PredictionModelsMonograph/files/")
    ## as <- get(load("../data/as.rda"))
    ## as[,age:=age5*5]
    ## pct1 <- prop.table(table(as$ct1))
    ## perg <- prop.table(table(as$erg.status))
    ## pgs <- prop.table(table(as$diaggs))
    ## l0 <- glm(erg.status~age5+diaggs+ct1+lmax+ppb5+lpsaden,data=as,family="binomial")
    ## as[asprogtime==0,asprogtime:=0.0001]
    ## G <- survival::survreg(Surv(asprogtime,asprog=="0")~1,data=as)
    ## F1 <- survival::survreg(Surv(asprogtime,asprog=="1")~erg.status+diaggs+ct1+lmax+ppb5+lpsaden,data=as)
    ## F2 <- survival::survreg(Surv(asprogtime,asprog==2)~erg.status+diaggs+ct1+lmax+ppb5+lpsaden,data=as)
    diaggs=diaggs.33=diaggs.34=ct1=erg.status=NULL
    m <- lava::lvm(~ age + lpsaden + ppb5 + lmax + ct1 + diaggs + erg.status)
    ## m <- categorical(m,~diaggs, K=3,p=pgs[-2])
    lava::distribution(m,~age) <- lava::normal.lvm(mean=65,sd=4)
    lava::distribution(m,~lmax) <- lava::normal.lvm(mean=3,sd=1.25)
    lava::distribution(m,~lpsaden) <- lava::normal.lvm(mean=-3,sd=0.9)
    lava::distribution(m,~ppb5) <- lava::normal.lvm(mean=2.9,sd=2)
    lava::distribution(m,~ct1) <- lava::binomial.lvm(p=0.092)
    lava::distribution(m,~erg.status) <- lava::binomial.lvm(p=.44)
    lava::distribution(m,~diaggs.33) <- lava::binomial.lvm(p=0.8)
    lava::distribution(m,~diaggs.34) <- lava::binomial.lvm(p=0.1)
    ## scale=- coef(G)[["(Intercept)"]]/G$scale
    ## shape = 1/G$scale
    lava::distribution(m,~t0) <- lava::coxWeibull.lvm(scale=exp(-1.69/.33),shape=1/.33)
    lava::distribution(m,~t1) <- lava::coxWeibull.lvm(scale=exp(-1.85/.9),shape=1/.9)
    lava::distribution(m,~t2) <- lava::coxWeibull.lvm(scale=exp(-1.6/.5),shape=1/.5)
    lava::regression(m,t1~erg.status+diaggs.33+diaggs.34+ct1+lmax+ppb5+lpsaden) <- -c(-0.9,0.6,-0.15,-0.32,-0.17,-0.1,-0.17)/0.9
    lava::regression(m,t2~erg.status+diaggs.33+diaggs.34+ct1+lmax+ppb5+lpsaden) <- -c(-0.04,0.18,0.66,0.14,0.04,-0.09,-0.03)/0.5
    m <- lava::eventTime(m,time~min(t0=0,t1=1,t2=2),"event")
    AS <- lava::sim(m,n)
    data.table::setDT(AS)
    AS[,diaggs:=factor(diaggs.33+diaggs.34,levels=c(0,1,2),labels=c("GNA","3/3","3/4"))]
    AS[,c("diaggs.33","diaggs.34","t0","t1","t2"):=NULL]
    AS[,ct1:=factor(ct1,levels=c("0","1"),labels=c("cT1","cT2"))]
    AS[,erg.status:=factor(erg.status,levels=c("0","1"),labels=c("neg","pos"))]
    AS[]
}


#----------------------------------------------------------------------
### simActiveSurveillance.R ends here
