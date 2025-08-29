library(lava)

set.seed(18)
m <- 500
n <- 500

###Binary Outcome
ciL <- numeric(m)
ciU <- numeric(m)
IPAest <- numeric(m)

for(i in 1:m){
  learndat <- sampleData(n,outcome="binary")
  testdat <- sampleData(n,outcome="binary")
  lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
  x=Score(list("LR(X1+X2+X7+X9)"=lr1),formula=Y~1,data=testdat, summary = "ipa")
  
  IPAest[i] <- x$IPA$score[[2,"IPA"]]
  ciL[i] <- x$IPA$score[[2,"lower"]]
  ciU[i] <- x$IPA$score[[2,"upper"]]
  print(i)
}
IPAest <- mean(IPAest)
coverage <- mean((IPAest>ciL)*(IPAest<ciU))
cat("Binary coverage = ", as.character(coverage))
### Binary outcome end

### Survival Outcome
ciL <- numeric(m)
ciU <- numeric(m)
IPAest <- numeric(m)

for(i in 1:m){
  trainSurv <- sampleData(n,outcome="survival")
  testdat <- sampleData(n,outcome="survival")
  cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
  x=Score(list("Cox(X1+X2+X7+X9)"=cox1),Surv(time,event)~1,data=testdat, summary = "ipa")
  
  IPAest[i] <- x$IPA$score[[2,"IPA"]]
  ciL[i] <- x$IPA$score[[2,"lower"]]
  ciU[i] <- x$IPA$score[[2,"upper"]]
  print(i)
}
IPAest <- mean(IPAest)
coverage <- mean((IPAest>ciL)*(IPAest<ciU))
cat("Survival coverage = ", as.character(coverage))
### Survival outcome end



ciL <- numeric(m)
ciU <- numeric(m)
IPAest <- numeric(m)

for(i in 1:m){
  ### Competing risks 
  trainCR <- sampleData(n,outcome="competing.risks")
  testCR <- sampleData(n,outcome="competing.risks")
  # Cause-specific Cox regression
  csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
  x = Score(list("CSC(X1+X2+X7+X9)"=csc1),formula=Hist(time,event)~1,data=testCR,summary = "ipa")
            
            
  
  IPAest[i] <- x$IPA$score[[2,"IPA"]]
  ciL[i] <- x$IPA$score[[2,"lower"]]
  ciU[i] <- x$IPA$score[[2,"upper"]]
  print(i)
}
IPAest <- mean(IPAest)
coverage <- mean((IPAest>ciL)*(IPAest<ciU))
cat("Competing risk coverage = ", as.character(coverage))










