rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 160000
W <- rnorm(n)
A_0 <- rexpit(W)
#Y <- rexpit(W + A)
L1_1 <- rbinom(n, size=1,prob=0.6)
L2_1 <- rbinom(n, size=1,prob=0.3)
A_1 <- rexpit(L1_1+L2_1+W)
Y1_1 <- rexpit(A_1+L1_1+L2_1+W)
Y2_1 <- rexpit(A_1+L1_1+L2_1+W) # teknisk set fejl
data_1 <- data.frame(W=W,A_0=A_0,L1_1=L1_1,L2_1=L2_1,A_1=A_1,Y1_1=Y1_1,Y2_1=Y2_1)
not_cens_1 <-Y1_1==0 & Y2_1 ==0

L1_2 <- rbinom(n, size=1,prob=0.1)
L2_2 <- rbinom(n, size=1,prob=0.7)
A_2 <- rexpit(L1_2+L2_2+W)
Y1_2 <- rexpit(A_2+L1_2+L2_2+W)
Y2_2 <- rexpit(A_2+L1_2+L2_2+W)
data_2_temp <- data.frame(L1_2=L1_2,L2_2=L2_2,A_2=A_2,Y1_2=Y1_2,Y2_2=Y2_2)
data_2_temp[!not_cens_1,] <- NA
not_cens_2 <- (Y1_2==0 & Y2_2 ==0) & not_cens_1
               
#data2 <- cbind(data_2_temp,data_1)
L1_3 <- rbinom(n, size=1,prob=0.7)
L2_3 <- rbinom(n, size=1,prob=0.7)
A_3 <- rexpit(L1_3+L2_3+W)
Y1_3 <- rexpit(A_3+L1_3+L2_3+W)
Y2_3 <- rexpit(A_3+L1_3+L2_3+W)
data_3_temp <- data.frame(L1_3=L1_3,L2_3=L2_3,A_3=A_3,Y1_3=Y1_3,Y2_3=Y2_3)
data_3_temp[!not_cens_2,] <- NA

data <- cbind(data_3_temp,data_2_temp,data_1)
#saveRDS(data,file"blabla.Rds")
library(riskRegression)
library(lava)
library(foreach)
d<-synthesizeLTMLE(data=data,A=c("A_0","A_1","A_2","A_3"),L=list(c("L1_1","L1_2","L1_3"),c("L2_1","L2_2","L2_3")),W="W",Y=list(c("Y1_1","Y1_2","Y1_3"),c("Y2_1","Y2_2","Y2_3")),time.points = 3)
d2<-sim(d,n=n)