library(ranger)
library(riskRegression)
library(survival)
data(GBSG2,package="pec")
fitranger <- ranger(Surv(time, cens)~horTh+age+pnodes+tsize+tgrade+estrec+progrec+menostat,
                    data = GBSG2,
                    oob.error=FALSE,
                    num.threads=1,
                    num.trees = 100)
scored <- Score(list(ranger=fitranger),data=GBSG2,times=1000,formula=Surv(time, cens)~1,conservative=TRUE,se.fit=FALSE,B=100,split.method="loob")
