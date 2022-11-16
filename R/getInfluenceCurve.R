getInfluenceCurve.AUC <- function(t,time,event, WTi, Wt, risk, MC, auc, nth.times, conservative, cens.model, one.step = FALSE, cv = FALSE){
  conservativeIFcalculation <- getIC0AUC(time,event,t,risk,WTi,Wt,auc)
  if (conservative[[1]] || cens.model[[1]] == "none"){
    conservativeIFcalculation[["ic0"]]
  }
  else {
    if (cens.model[[1]] == "KaplanMeier"){
      conservativeIFcalculation[["ic0"]]+getInfluenceFunctionAUCKMCensoringTerm(time,event,t,conservativeIFcalculation[["ic0Case"]],
                                                                                conservativeIFcalculation[["ic0Control"]], conservativeIFcalculation[["weights"]],
                                                                                conservativeIFcalculation[["firsthit"]],conservativeIFcalculation[["muCase"]],
                                                                                conservativeIFcalculation[["muControls"]], conservativeIFcalculation[["nu"]], Wt[1], auc)
    }
    else if (cens.model[[1]] == "cox" || cens.model[[1]] == "discrete"){
      n <- length(time)
      ic0CaseOld <- conservativeIFcalculation[["ic0Case"]]
      ic0ControlOld <- conservativeIFcalculation[["ic0Control"]]
      Phi <- conservativeIFcalculation[["muCase"]] * conservativeIFcalculation[["muControls"]] / (n*n)
      weights <- conservativeIFcalculation[["weights"]]
      aucLPO <- auc
      cases.index <- time<=t & event==1
      controls.index1 <- time>t
      controls.index2 <-  event==2 & time <= t
      controls.index <- controls.index1 | controls.index2
      w.cases <- weights[cases.index]
      w.controls <- weights[controls.index]
      ic0Case <- ic0CaseOld[cases.index]*w.cases
      ic0Control <- ic0ControlOld[controls.index]*w.controls
      if (!MC$saveCoxMemory){
        ic.weights <- MC[[2]][[nth.times]] ## load IF from Censoring weights
        ## Part of influence function related to Weights
        ic.weightsCase <- as.numeric(rowSumsCrossprod(as.matrix(ic0Case), ic.weights[cases.index,], 0)) ## part of second term
        ic.weightsControl <- as.numeric(rowSumsCrossprod(as.matrix(ic0Control), ic.weights[controls.index,], 0)) ## part of third term
        ic.weightsCC <- (1/(Phi*n^2))*(ic.weightsCase+ic.weightsControl)
        ## Note that the first and fourth term in IF_mu*nu(Q)/mu(Q)^2 = AUC_tau / mu(Q) * IF_mu is IF.AUC0 in the below
        ## ## Part of influence function related to Phi, i.e. mu(Q) 
        icPhiCase <- as.numeric(rowSumsCrossprod(as.matrix(w.cases),ic.weights[cases.index,],0)) ## part of second term
        icPhiControl <- as.numeric(rowSumsCrossprod(as.matrix(w.controls),ic.weights[controls.index,],0))## part of third term
        icPhi <- (aucLPO/Phi)*(((1/n)*icPhiCase)*(1/n)*sum(w.controls)+((1/n)*icPhiControl)*(1/n)*sum(w.cases))
      }
      else {
        wdata <- MC[[3]]
        fit <- MC[[2]]
        TiMinus <- MC[[4]]
        ic0beforet <- ic0CaseOld*weights*cases.index+ic0ControlOld*weights*controls.index2
        ic0aftert <- ic0ControlOld*weights*controls.index1
        ic.weightsCC <- (1/(Phi*n^2))*(predictCoxWeights(fit, diag=TRUE,newdata = wdata, times = TiMinus,weights=ic0beforet)+predictCoxWeights(fit, diag=FALSE,newdata = wdata,times = t,weights=ic0aftert))
        weightsbeforet <- (cases.index)*weights*(1/n)*sum(w.controls)+ (controls.index2)*weights*(1/n)*sum(w.cases)
        weightsaftert <- (controls.index1)*weights*(1/n)*sum(w.cases)
        icPhiBeforet <- predictCoxWeights(fit, diag=TRUE,newdata = wdata, times = TiMinus,weights=weightsbeforet)
        icPhiAftert <- predictCoxWeights(fit, diag=FALSE,newdata = wdata, times = t,weights=weightsaftert)
        icPhi <- (aucLPO/Phi)*((icPhiBeforet+icPhiAftert)*(1/n))
      }
      conservativeIFcalculation[["ic0"]]+ic.weightsCC-icPhi
    }
    else {
      stop("Censoring model not yet implemented. ")
    }
  }
}

getInfluenceCurve.Brier <- function(t,
                                    time,
                                    IC0,
                                    residuals,
                                    WTi,
                                    Wt,
                                    IC.G,
                                    cens.model,
                                    nth.times=NULL,
                                    conservative,
                                    event,
                                    one.step = FALSE){
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
    if (conservative[[1]] || cens.model[[1]] == "none"){
      IC0
    }
    else if (cens.model[[1]]=="cox" || cens.model[[1]] == "discrete") {## Cox regression
        if (!IC.G$saveCoxMemory){
          ic.weights <- IC.G[[2]][[nth.times]]
          IF.Brier <- IC0+as.numeric(rowSumsCrossprod(as.matrix(residuals),ic.weights,0)) / N
        }
        else {
          wdata <- IC.G[[3]]
          fit <- IC.G[[2]]
          TiMinus <- IC.G[[4]]
          res1 <- (time <= t & event != 0)*residuals
          res2 <- (time > t)*residuals
          IF.Brier <- IC0 + (predictCoxWeights(fit,newdata = wdata,times = TiMinus,diag=TRUE,weights=res1) + predictCoxWeights(fit,newdata = wdata,times = t,diag=FALSE,weights=res2))/N
        }
        IF.Brier
    }else if (cens.model[[1]] == "KaplanMeier"){
        IC0 + getInfluenceFunctionBrierKMCensoringTerm(t,time,residuals,event)
    }
    else {
      stop("Non conservative options with cens.model not being a Cox model or KaplanMeier are not implemented. ")
    }
}

getInfluenceCurve.Brier.covariates <- function(tau,time,residuals,risk,status,GTiminus,Gtau,IC.data) {
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
    IC[i] <- IC.C.term
  }
  IC
}

getInfluenceCurve.AUC.covariates <- function(t,n,time,status,risk,GTiminus,Gtau,AUC,IC.data){ ## needs to be fixed later
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
