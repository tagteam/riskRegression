getInfluenceCurveHelper <- function(time,status,tau,risk,GTiminus,Gtau,AUC){
  urisk <- unique(risk)
  if (length(urisk)==length(risk)){
    getInfluenceFunctionAUCKMCensoring(time,status,tau,risk,GTiminus,Gtau,AUC,FALSE)
  }
  else {
    n <- length(time)
    nutauParti.ties <-rep(NA,n)
    for (val in urisk){
      indexes <- which(risk==val)
      nutauParti.ties[indexes] <- mean(1*(risk == val & time <= tau & status == 1)/GTiminus)-1*(risk[indexes] == val & time[indexes]<= tau & status[indexes] == 1)/(n*GTiminus[indexes])
    }
    numAUCties <- mean(nutauParti.ties * 1*(time > tau)/Gtau) + mean(nutauParti.ties * 1*(time <= tau & status == 2)/GTiminus)
    denAUC <- mean(1*(time <= tau & status == 1)/GTiminus)*mean(time > tau)/Gtau + mean(1*(time <= tau & status == 1)/GTiminus)*mean(1*(time <= tau & status == 2)/GTiminus)
    AUC.ties.part<-(numAUCties)/denAUC
    AUC.noties <- AUC-0.5*AUC.ties.part
    IF.noties <- getInfluenceFunctionAUCKMCensoring(time,status,tau,risk,GTiminus,Gtau,AUC.noties,FALSE)
    IF.ties <- getInfluenceFunctionAUCKMCensoring(time,status,tau,risk,GTiminus,Gtau,AUC.ties.part,TRUE)
    IF.noties+0.5*IF.ties
  }
}

getInfluenceCurve.AUC.cox <- function(t,time,event, WTi, Wt, risk, ID, MC, nth.times){
  n <- length(time)
  cases.index <- time<=t & event==1
  controls.index1 <- time>t
  controls.index2 <-  event==2 & time <= t
  controls.index <- controls.index1 | controls.index2
  
  cc.status <- factor(rep("censored",n),levels=c("censored","case","control"))
  cc.status[cases.index] <- "case"
  cc.status[controls.index] <- "control"
  
  id.cases <- ID[cc.status=="case"]
  id.controls <- ID[cc.status=="control"]
  id.censored <- ID[cc.status=="censored"]
  
  weights.cases <- cases.index/WTi
  weights.controls <- 1/Wt * controls.index1  + 1/WTi * controls.index2
  w.cases <-weights.cases[cases.index]
  w.controls <- weights.controls[controls.index]
  weightMatrix <- outer(w.cases, w.controls, "*")
  Phi <- (1/n^2)*sum(w.cases)*sum(w.controls)
  n.cases <- sum(cases.index)
  n.controls <- sum(controls.index)
  risk.cases <- risk[cases.index]
  risk.controls <- risk[controls.index]
  ic0Case <- rep(0,n.cases)
  ic0Control <- rep(0,n.controls)
  for (i in 1:n.controls){
    ic0Case <- ic0Case + (1*(risk.cases > risk.controls[i])+0.5*1*(risk.cases == risk.controls[i]))*w.cases*w.controls[i]
  }
  for (i in 1:n.cases){
    ic0Control <- ic0Control + (1*(risk.cases[i] > risk.controls)+0.5*1*(risk.cases[i] == risk.controls))*w.cases[i]*w.controls
  }
  # auc <- AUCijFun(risk.cases,risk.controls)
  # auc <- auc*weightMatrix
  # auc[is.na(auc)] <- 0
  nu1tauPm <- (1/n^2)*sum(ic0Case)
  aucLPO <- nu1tauPm*(1/Phi)
  # ic0Case <- rowSums(auc)
  # ic0Control <- colSums(auc)
  ic0 <- (1/(Phi*n))*c(ic0Case, ic0Control)-2*aucLPO
  aucDT <- data.table(ID = c(id.cases,id.controls,id.censored), IF.AUC0 = c(ic0, rep(-2*aucLPO,length(id.censored))))
  ic.weights <- matrix(0,n,n)
  k=0 ## counts subject-times with event before t
  for (i in 1:n){
    if (i %in% id.cases){
      k=k+1
      ic.weights[i,] <- MC$IC.subject[i,k,]/WTi[i]
    }
    else if (i %in% id.controls){ 
      if (time[i] > t){## min(T,C)>t
        ic.weights[i,] <- MC$IC.times[i,nth.times,]/Wt[i]
      }
      else { ## min(T,C)<= t and status == 2
        k=k+1
        ic.weights[i,] <- MC$IC.subject[i,k,]/WTi[i]
      }
    }
  }
  ## ## Part of influence function related to Weights
  ic.weightsCase <- as.numeric(rowSumsCrossprod(as.matrix(ic0Case), ic.weights[cases.index,], 0))
  ic.weightsControl <- as.numeric(rowSumsCrossprod(as.matrix(ic0Control), ic.weights[controls.index,], 0))
  ic.weightsCC <- (1/(Phi*n^2))*(ic.weightsCase+ic.weightsControl)
  ## ## Part of influence function related to Phi
  ## icPhiCase <- colMeans(ic.weights[which.cases,])
  icPhiCase <- as.numeric(rowSumsCrossprod(as.matrix(w.cases),ic.weights[cases.index,],0))
  icPhiControl <- as.numeric(rowSumsCrossprod(as.matrix(w.controls),ic.weights[controls.index,],0))
  icPhi <- (aucLPO/Phi)*((weights.cases-(1/n)*icPhiCase)*(1/n)*sum(weights.controls)+(weights.controls-(1/n)*icPhiControl)*(1/n)*sum(weights.cases)) - 2*aucLPO
  ## ## Combine all parts of influence function
  ## ic1 <- data.table(ID=data[["ID"]], "ic.weightsCC" = ic.weightsCC, "icPhi" = icPhi)
  data.table::setkey(aucDT,ID)
  aucDT[["IF.AUC0"]]-ic.weightsCC-icPhi
}

getInfluenceCurve.AUC.conservative <- function(t,time,event, WTi, Wt, risk, ID){
  n <- length(time)
  cases.index <- time<=t & event==1
  controls.index1 <- time>t
  controls.index2 <-  event==2 & time <= t
  controls.index <- controls.index1 | controls.index2
  
  cc.status <- factor(rep("censored",n),levels=c("censored","case","control"))
  cc.status[cases.index] <- "case"
  cc.status[controls.index] <- "control"
  
  id.cases <- ID[cc.status=="case"]
  id.controls <- ID[cc.status=="control"]
  id.censored <- ID[cc.status=="censored"]
  
  weights.cases <- cases.index/WTi
  weights.controls <- 1/Wt * controls.index1  + 1/WTi * controls.index2
  w.cases <-weights.cases[cases.index]
  w.controls <- weights.controls[controls.index]
  weightMatrix <- outer(w.cases, w.controls, "*")
  Phi <- (1/n^2)*sum(w.cases)*sum(w.controls)
  n.cases <- sum(cases.index)
  n.controls <- sum(controls.index)
  risk.cases <- risk[cases.index]
  risk.controls <- risk[controls.index]
  ic0Case <- rep(0,n.cases)
  ic0Control <- rep(0,n.controls)
  for (i in 1:n.controls){
    ic0Case <- ic0Case + (1*(risk.cases > risk.controls[i])+0.5*1*(risk.cases == risk.controls[i]))*w.cases*w.controls[i]
  }
  for (i in 1:n.cases){
    ic0Control <- ic0Control + (1*(risk.cases[i] > risk.controls)+0.5*1*(risk.cases[i] == risk.controls))*w.cases[i]*w.controls
  }
  # auc <- AUCijFun(risk.cases,risk.controls)
  # auc <- auc*weightMatrix
  # auc[is.na(auc)] <- 0
  nu1tauPm <- (1/n^2)*sum(ic0Case)
  aucLPO <- nu1tauPm*(1/Phi)
  # ic0Case <- rowSums(auc)
  # ic0Control <- colSums(auc)
  ic0 <- (1/(Phi*n))*c(ic0Case, ic0Control)-2*aucLPO
  aucDT <- data.table(ID = c(id.cases,id.controls,id.censored), IF.AUC0 = c(ic0, rep(-2*aucLPO,length(id.censored))))
  ## ## Part of influence function related to Phi
  ## icPhiCase <- colMeans(ic.weights[which.cases,])
  icPhi <- (aucLPO/Phi)*((weights.cases)*(1/n)*sum(weights.controls)+(weights.controls)*(1/n)*sum(weights.cases)) - 2*aucLPO
  ## ## Combine all parts of influence function
  ## ic1 <- data.table(ID=data[["ID"]], "ic.weightsCC" = ic.weightsCC, "icPhi" = icPhi)
  data.table::setkey(aucDT,ID)
  aucDT[["IF.AUC0"]]-icPhi
}

getInfluenceCurve.Brier <- function(t,
                                    time,
                                    IC0,
                                    residuals,
                                    WTi,
                                    Wt,
                                    IC.G,
                                    cens.model,
                                    nth.times=NULL){
    ##
    ## Compute influence function of Brier score estimator using weights of the reverse Cox model
    ## This function evaluates the part of influence function which is related to the IPCW weights
    ## The other part is IC0.
    ##
    ## \frac{1}{n}\sum_{i=1}^n
    ## m_{t,n}^{(1)}(X_i)
    ## [\frac{I_{T_i\leq t}\Delta_i}{G^2(T_i\vert Z_i)}IF_G(T_i,X_k; X_i)+\frac{I_{T_i>t}}{G^2(t|Z_i)}IF_G(t,X_k; X_i)]
    ## with
    ## IF_G(t,X_k; X_i)=-\exp(-\Lambda(t\vert Z_i))IF_{\Lambda}(t,X_k; X_i)
    ##
    ## IC_G(t,z;x_k) is an array with dimension (nlearn=N, gtimes, newdata)
    ## where gtimes = subject.times (Weights$IC$IC.subject) or times (Weights$IC$IC.times)
    ## and subject.times=Y[(((Y<=max(times))*status)==1)]
    ##
    ## don't square the weights because they will be multiplied with the
    ## residuals that are already weighted
    ##
    N <- length(residuals)
    if (cens.model=="cox") {## Cox regression
        ic.weights <- matrix(0,N,N)
        k=0 ## counts subject-times with event before t
        for (i in 1:N){
            if (residuals[i]>0){
                if (time[i]<=t){ ## min(T,C)<=t, note that (residuals==0) => (status==0)
                    k=k+1
                    ic.weights[i,] <- IC.G$IC.subject[i,k,]/(WTi[i])
                }else{## min(T,C)>t
                    ic.weights[i,] <- IC.G$IC.times[i,nth.times,]/(Wt[i])
                }
            }
        }
        IF.Brier <- ic.weights*residuals
        IF.Brier <- IC0-colMeans(IF.Brier)
        IF.Brier
    }else{
        stop("Use C++ function instead. ")
    }
}

#NOTE: Requires cox model for censoring (and time)!
getInfluenceCurve.Brier.covariates <- function(tau,time,risk,status,GTiminus,Brier,IC.data) {
  n <- length(time)
  IC <- rep(NA,n)
  fit.time <- IC.data$fit.time
  fit.cens <- IC.data$fit.cens
  wdata<-IC.data$wdata
  ## term involving f_i(t,z) is 
  ## $$
  ## (1-2R(\tau |Z_i))\left(\frac{I(\tilde{T}_i \leq \tau, \Delta_i = 0)}{G(\tilde{T}_i|Z_i)S(\tilde{T}_i|Z_i)}(S(\tilde{T}_i|Z_i)-S(\tau|Z_i))-\int_0^{\tilde{T}_i \wedge \tau} \frac{(S(s|Z_i)-S(\tau|Z_i))}{G(s|Z_i)^2S(s|Z_i)^2}P(ds,0|Z_i)\right)
  ## $$
  predCens <- predictCox(fit.cens, time,newdata = wdata)
  Gtimes <- predCens$survival
  cumhazardCXi <- predCens$cumhazard
  Stimes <- predictCox(fit.time, time,newdata = wdata)$survival
  is.comprisk <- !is.null(IC.data$fitCSC)
  if (is.comprisk){
    fitCSC <- IC.data$fitCSC
    F1 <- predictRisk(fitCSC,newdata=wdata,times=time,cause=1)
    F1tau <- c(predictRisk(fitCSC,newdata=wdata,times=tau,cause=1))
  }
  else {
    Stau <- rms::survest(fit.time,newdata=wdata,times=tau,se.fit=FALSE)$surv
    Stau <- unname(Stau)
  }
  for (i in 1:n){
    cum <- cumhazardCXi[i,]
    jumps <- diff(c(0,cum))
    if (!is.comprisk){
      IC.C.term <- (1-2*risk[i])*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(Stimes[i,i]-Stau[i])-sum( 1*(time <= tau & time <= time[i])*((Stimes[i,]-Stau[i])*jumps / (Gtimes[i,]*Stimes[i,]))))
    }
    else {
      IC.C.term <- (1-2*risk[i])*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F1tau[i]-F1[i,i])-sum( 1*(time <= tau & time <= time[i])*((F1tau[i]-F1[i,])*jumps / (Gtimes[i,]*Stimes[i,]))))
    }
    # IC.C.term <- 0
    IC[i] <- 1*(time[i] <= tau & status[i] == 1 )* (1-2*risk[i]) * 1/GTiminus[i] + IC.C.term + risk[i]^2 - Brier
  }
  IC
}

getInfluenceCurve.Brier.covariates.use.squared <- function(tau,time,residuals,risk,status,GTiminus,Gtau,IC.data) {
  n <- length(time)
  IC <- rep(NA,n)
  fit.time <- IC.data$fit.time
  fit.cens <- IC.data$fit.cens
  wdata<-IC.data$wdata
  Brier <- mean(residuals)
  predCens <- predictCox(fit.cens, time,newdata = wdata)
  Gtimes <- predCens$survival
  cumhazardCXi <- predCens$cumhazard
  Stimes <- predictCox(fit.time, time,newdata = wdata)$survival
  is.comprisk <- !is.null(IC.data$fitCSC)
  Stau <- rms::survest(fit.time,newdata=wdata,times=tau,se.fit=FALSE)$surv
  Stau <- unname(Stau)
  if (is.comprisk){
    fitCSC <- IC.data$fitCSC
    F1 <- predictRisk(fitCSC,newdata=wdata,times=time,cause=1)
    F1tau <- c(predictRisk(fitCSC,newdata=wdata,times=tau,cause=1))
  }
  for (i in 1:n){
    cum <- cumhazardCXi[i,]
    jumps <- diff(c(0,cum))
    IC.C.term.part <- risk[i]^2 * Stau[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])-sum( 1*(time <= tau & time <= time[i])*jumps / (Gtimes[i,]*Stimes[i,])))
    if (!is.comprisk){
      IC.C.term <- IC.C.term.part+ (1-risk[i])^2*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(Stimes[i,i]-Stau[i])-sum( 1*(time <= tau & time <= time[i])*((Stimes[i,]-Stau[i])*jumps / (Gtimes[i,]*Stimes[i,]))))
    }
    else {
      IC.C.term <- IC.C.term.part + (1-risk[i])^2*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F1tau[i]-F1[i,i])-sum( 1*(time <= tau & time <= time[i])*((F1tau[i]-F1[i,])*jumps / (Gtimes[i,]*Stimes[i,]))))
    }
    # IC.C.term <- 0
    IC[i] <- residuals[i]-Brier + IC.C.term
  }
  IC
}

getInfluenceCurve.AUC.covariates <- function(t,n,time,status,risk,GTiminus,Gtau,AUC,IC.data){
    tau <- t
    X <- risk
    #estimate int 1{X_i > x, t' > tau} dP(t',x)/G(tau | x')
    int1nu <- rep(NA,n)
    #estimate int 1{t' > tau} dP(t',x)/G(tau | x')
    int1mu <- mean(1*(time > tau)/(Gtau))
    #estimate int 1{X_i < x, t <= tau} dP(t,1,x)/G(t | x)
    int2nu <- rep(NA,n)
    #estimate int 1{t <= tau} dP(t,1,x)/G(t | x)
    int2mu <- mean(1*(time <= tau & status == 1)/(GTiminus))
    #estimate int 1{X_i > x, t <= tau} dP(t,1,x)/G(t | x)
    int3nu <- rep(NA,n)
    #estimate int 1{t <= tau} dP(t,1,x)/G(t | x)
    int3mu <-mean(1*(time <= tau & status == 2)/(GTiminus))

    for (i in 1:n){
        int1nu[i] <- mean(1*(X[i] > X & time > tau)/(Gtau))
        int2nu[i] <- mean(1*(X[i] < X & time <= tau & status == 1)/(GTiminus))
        int3nu[i] <- mean(1*(X[i] > X & time <= tau & status == 2)/(GTiminus))
    }

    ic <- rep(0,n)
    if (any(Gtau == 0) || any(GTiminus == 0)) {
        stop("Some censoring weights are 0. Pick another censoring model or retry with a larger data set.")
    }
    #main loop
    fhat.tau <- rep(0,n)
    fhat.Ti <- rep(0,n)
    mu1 <- int2mu * int1mu + int2mu*int3mu
    nu1 <- AUC*mu1
    
    fit.time <- IC.data$fit.time
    fit.cens <- IC.data$fit.cens
    wdata<-IC.data$wdata
    ## term involving f_i(t,z) is 
    ## $$
    ## (1-2R(\tau |Z_i))\left(\frac{I(\tilde{T}_i \leq \tau, \Delta_i = 0)}{G(\tilde{T}_i|Z_i)S(\tilde{T}_i|Z_i)}(S(\tilde{T}_i|Z_i)-S(\tau|Z_i))-\int_0^{\tilde{T}_i \wedge \tau} \frac{(S(s|Z_i)-S(\tau|Z_i))}{G(s|Z_i)^2S(s|Z_i)^2}P(ds,0|Z_i)\right)
    ## $$
    predCens <- predictCox(fit.cens, time,newdata = wdata)
    Gtimes <- predCens$survival
    cumhazardCXi <- predCens$cumhazard
    Stimes <- predictCox(fit.time, time,newdata = wdata)$survival
    is.comprisk <- !is.null(IC.data$fitCSC)
    Stau <- unname(rms::survest(fit.time,newdata=wdata,times=tau,se.fit=FALSE)$surv)
    if (is.comprisk){
      fitCSC <- IC.data$fitCSC
      F1 <- predictRisk(fitCSC,newdata=wdata,times=time,cause=1)
      F1tau <- c(predictRisk(fitCSC,newdata=wdata,times=tau,cause=1))
      F2 <- predictRisk(fitCSC,newdata=wdata,times=time,cause=2)
      F2tau <- c(predictRisk(fitCSC,newdata=wdata,times=tau,cause=2))
    }

    for (i in 1:n){
        cum <- cumhazardCXi[i,]
        jumps <- diff(c(0,cum))
        
        # if (any(fhat.Ti != 0) || any(fhat.tau != 0)) browser()
        # #calculate fhat(tau,X_i) for i = 1, ..., n
        term1nu <- 1*(time[i] <= tau & status[i] == 1)/GTiminus[i] * int1nu[i]
        if (!is.comprisk){
          term2nu <- int1nu[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(Stimes[i,i]-Stau[i])-sum( 1*(time <= tau & time <= time[i])*((Stimes[i,]-Stau[i])*jumps / (Gtimes[i,]*Stimes[i,])))) -mean(int1nu * 1*(time <= tau & status == 1) * 1/GTiminus)
          term6nu <- 0
          term8nu <- int3nu[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(Stimes[i,i]-Stau[i])-sum( 1*(time <= tau & time <= time[i])*((Stimes[i,]-Stau[i])*jumps / (Gtimes[i,]*Stimes[i,])))) -mean(int3nu * 1*(time <= tau & status == 1) * 1/GTiminus)
          intmu <- 1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(Stimes[i,i]-Stau[i])-sum( 1*(time <= tau & time <= time[i])*((Stimes[i,]-Stau[i])*jumps / (Gtimes[i,]*Stimes[i,])))- mean(1*(time <= tau & status == 1) * 1/GTiminus)
          term6mu <- 0
        }
        else {
          term2nu <- int1nu[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F1tau[i]-F1[i,i])-sum( 1*(time <= tau & time <= time[i])*((F1tau[i]-F1[i,])*jumps / (Gtimes[i,]*Stimes[i,])))) -mean(int1nu * 1*(time <= tau & status == 1) * 1/GTiminus)
          term6nu <- int2nu[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F2tau[i]-F2[i,i])-sum( 1*(time <= tau & time <= time[i])*((F2tau[i]-F2[i,])*jumps / (Gtimes[i,]*Stimes[i,])))) -mean(int2nu * 1*(time <= tau & status == 2) * 1/GTiminus)
          term8nu <- int3nu[i] * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F1tau[i]-F1[i,i])-sum( 1*(time <= tau & time <= time[i])*((F1tau[i]-F1[i,])*jumps / (Gtimes[i,]*Stimes[i,])))) -mean(int3nu * 1*(time <= tau & status == 1) * 1/GTiminus)
          intmu <- 1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F1tau[i]-F1[i,i])-sum( 1*(time <= tau & time <= time[i])*((F1tau[i]-F1[i,])*jumps / (Gtimes[i,]*Stimes[i,])))- mean(1*(time <= tau & status == 1) * 1/GTiminus)
          term6mu <- int2mu * (1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])*(F2tau[i]-F2[i,i])-sum( 1*(time <= tau & time <= time[i])*((F2tau[i]-F2[i,])*jumps / (Gtimes[i,]*Stimes[i,])))- mean(1*(time <= tau & status == 2) * 1/GTiminus))
        }
        term3nu <- 1*(time[i] > tau)/Gtau[i] * int2nu[i]
        term4nu <- int2nu[i] * Stau[i]*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])-sum( 1*(time <= tau & time <= time[i])*jumps / (Gtimes[i,]*Stimes[i,]))) - mean(int2nu * 1*(time > tau) * 1/Gtau)
        term5nu <- 1*(time[i] <= tau & status[i] == 2)/GTiminus[i] * int2nu[i]
        term7nu <- 1*(time[i] <= tau & status[i] == 1)/GTiminus[i] * int3nu[i]
        IFnu <- term1nu+term2nu+term3nu+term4nu+term5nu+term6nu+term7nu+term8nu

        term1mu <- 1*(time[i] <= tau & status[i] == 1)/GTiminus[i] * int1mu
        term2mu <- int1mu * intmu
        term3mu <- 1*(time[i] > tau)/Gtau[i] * int2mu
        term4mu <- int2mu*(Stau[i]*(1*(status[i] == 0 & time[i] <= tau)/(Gtimes[i,i]*Stimes[i,i])-sum( 1*(time <= tau & time <= time[i])*jumps / (Gtimes[i,]*Stimes[i,]))) - mean(1*(time > tau) * 1/Gtau))
        term5mu <- 1*(time[i] <= tau & status[i] == 2)/GTiminus[i] * int2mu
        term7mu <- 1*(time[i] <= tau & status[i] == 1)/GTiminus[i] * int3mu
        term8mu <- int3mu*intmu
        IFmu <- term1mu+term2mu+term3mu+term4mu+term5mu+term6mu+term7mu+term8mu
        ic[i] <- (IFnu * mu1 - nu1 * IFmu)/(mu1^2)
    }
    ic
}
