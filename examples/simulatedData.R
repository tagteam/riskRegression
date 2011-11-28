library(riskRegression)
library(prodlim)
d <- prodlim:::SimCompRisk(100)
check.code("prodlim")
check.code("riskRegression")

ttt <- seq(0,80,1)
f <- riskRegression(Hist(time,cause)~X1+X2,times=ttt,link="relative",data=d)
b <- ARR(Hist(time,cause)~X1+X2,times=ttt,data=d)

b2 <- ARR(Hist(time,cause)~pow(X1)+X2,times=ttt,data=d)

fg <- riskRegression(Hist(time,cause)~X1+X2,times=ttt,link="prop",data=d)
fl <- LRR(Hist(time,cause)~X1+X2,times=ttt,data=d)

## set.seed(17)
ttt <- seq(0,4,.1)
nix <- do.call("rbind",lapply(1:10,function(x){
  print(x)
  d <-  SimCompRisk(1000,cova.X=list("rbinom",size=1,prob=.7),cuminc1.coef=2,cuminc2.coef=1)
  ## table(d$cause,d$X)
  f <- riskRegression(Hist(time,cause)~X,times=ttt,data=d)
  fg <- riskRegression(Hist(time,cause)~X,times=ttt,data=d,link="prop")
  cbind(exp(f$timeConstCoef),exp(fg$timeConstCoef))
}))
nix
colMeans(nix)

## set.seed(17)
d <-  SimCompRisk(10000,cova.X=list("rbinom",size=1,prob=.7),cuminc1.coef=2,cuminc2.coef=1)
table(d$cause,d$X)
u <- prodlim(Hist(time,cause)~X,data=d)
plot(u,xlim=c(0,2))
ppp <- predict(u,newdata=data.frame(X=0:1),times=2)
ppp[[2]]/ppp[[1]]
