getInfluenceCurveHelper <- function(time,status,tau,risk,GTiminus,Gtau,AUC){
  urisk <- unique(risk)
  if (length(urisk)==length(risk)){
    getInfluenceFunctionAUC(time,status,tau,risk,GTiminus,Gtau,AUC,FALSE,FALSE,FALSE)
  }
  else {
    n <- length(time)
    # nutauParti<-rep(NA,n)
    nutauParti.ties <-rep(NA,n)
    for (val in urisk){
      indexes <- which(risk==val)
      nutauParti.ties[indexes] <- mean(1*(risk == val & time <= tau & status == 1)/GTiminus)-1*(risk[indexes] == val & time[indexes]<= tau & status[indexes] == 1)/(n*GTiminus[indexes])
    }
    # for (i in 1:n){
    #     nutauParti[i] <- mean(1*(risk > risk[i] & time <= tau & status == 1)/GTiminus)
    #     riskTies <- risk
    #     riskTies[i] <- -1
    #     nutauParti.ties[i] <- mean(1*(riskTies == risk[i] & time <= tau & status == 1)/GTiminus)
    # }
    # numAUC <- mean(nutauParti * 1*(time > tau)/Gtau) + mean(nutauParti * 1*(time <= tau & status == 2)/GTiminus)
    numAUCties <- mean(nutauParti.ties * 1*(time > tau)/Gtau) + mean(nutauParti.ties * 1*(time <= tau & status == 2)/GTiminus)
    denAUC <- mean(1*(time <= tau & status == 1)/GTiminus)*mean(time > tau)/Gtau + mean(1*(time <= tau & status == 1)/GTiminus)*mean(1*(time <= tau & status == 2)/GTiminus)
    AUC.ties.part<-(numAUCties)/denAUC
    AUC.noties <- AUC-0.5*AUC.ties.part
    IF.noties <- getInfluenceFunctionAUC(time,status,tau,risk,GTiminus,Gtau,AUC.noties,FALSE,FALSE,FALSE)
    IF.ties <- getInfluenceFunctionAUC(time,status,tau,risk,GTiminus,Gtau,AUC.ties.part,FALSE,TRUE,FALSE)
    IF.noties+0.5*IF.ties
  }
}
