getInfluenceFunction.AUC.censoring.term <- function(time,event,t, IFcalculationList, MC, cens.model, Wt, auc,nth.times){
  if (cens.model[[1]] == "KaplanMeier"){
    ind.controls<-rep(NA,length(time))
    controls.index1 <- IFcalculationList[["controls1"]]
    controls.index2 <- IFcalculationList[["controls2"]]
    ind.controls[controls.index1] <- 1
    ind.controls[controls.index2] <- 0
    start.controls1 <- sindex(ind.controls[controls.index1 | controls.index2],0)
    getInfluenceFunctionAUCKMCensoringTerm(time,
                                           event,
                                           t,
                                           IFcalculationList[["ic0Case"]],
                                           IFcalculationList[["ic0Control"]],
                                           IFcalculationList[["weights"]],
                                           IFcalculationList[["firsthit"]],
                                           IFcalculationList[["muCase"]],
                                           IFcalculationList[["muControls"]], 
                                           IFcalculationList[["nu"]],
                                           Wt[1],
                                           auc,
                                           start.controls1)
  }
  else if (cens.model[[1]] == "cox"){
    n <- length(time)
    ic0Case <- IFcalculationList[["ic0Case"]]
    ic0Control <- IFcalculationList[["ic0Control"]]
    Phi <- IFcalculationList[["muCase"]] * IFcalculationList[["muControls"]] / (n*n)
    weights <- IFcalculationList[["weights"]]
    muCase <- IFcalculationList[["muCase"]]
    muControls <- IFcalculationList[["muControls"]]
    aucLPO <- auc
    cases.index <- IFcalculationList[["cases"]]
    controls.index <- IFcalculationList[["controls"]]
    w.cases <- weights[cases.index]
    w.controls <- weights[controls.index]
    if (!MC$censoring.save.memory){
      ic.weights <- MC[[2]][[nth.times]] ## load IF from Censoring weights
      icPart <- as.numeric(rowSumsCrossprod(as.matrix(1/(Phi*n^2)*ic0Case-(aucLPO/Phi)*(1/n^2)*muControls*w.cases), ic.weights[cases.index,], 0)) +
        as.numeric(rowSumsCrossprod(as.matrix(1/(Phi*n^2)*ic0Control-(aucLPO/Phi)*(1/n^2)*muCase*w.controls),ic.weights[controls.index,],0))
    }
    else {
      wdata <- MC[[3]]
      fit <- MC[[2]]
      TiMinus <- MC[[4]]
      ic0CaseOld <- rep(0,n)
      ic0CaseOld[cases.index] <- ic0Case
      ic0ControlOld <- rep(0,n)
      ic0ControlOld[controls.index] <- ic0Control
      controls.index1 <- IFcalculationList[["controls1"]]
      controls.index2 <- IFcalculationList[["controls2"]]
      Wbeforet <- (1/(Phi*n^2))*(ic0CaseOld*cases.index+ic0ControlOld*controls.index2)-
        (1/n)*(aucLPO/Phi)*((cases.index)*weights*(1/n)*muControls + (controls.index2)*weights*(1/n)*muCase)
      Waftert <- (1/(Phi*n^2))*ic0ControlOld*controls.index1- 
        (1/n)*(aucLPO/Phi)*(controls.index1)*weights*(1/n)*muCase
      ## First term gives for i'th entry: 1/n \sum_j weights[j] * \hat{f}_i(\tilde{T}_j-,X_j); 
      ## Next one does: 1/n \sum_j weights[j] * \hat{f}_i(tau,X_j) for Cox
      icPart <- predictCoxWeights(fit, diag=TRUE,newdata = wdata, times = TiMinus,weights=Wbeforet, isBeforeTau = TRUE, tau = t)+
        predictCoxWeights(fit, diag=FALSE,newdata = wdata,times = t,weights=Waftert)
    }
    icPart
  }
  else {
    warning("Censoring model not yet implemented. Reverting to conservative = TRUE for AUC. ")
    0
  }
}

getInfluenceCurve.AUC <- function(t,time,event, WTi, Wt, risk, MC, auc, nth.times, conservative, cens.model){
  ## assert that time is sorted in ascending order with stop
  if (is.unsorted(time)){
    stop("Internal error. Time is not sorted in ascending order. ")
  }
  
  conservativeIFcalculation <- getIC0AUC(time,event,t,risk,WTi,Wt,auc)
  if (conservative[[1]] || cens.model[[1]] == "none"){
    conservativeIFcalculation[["ic0"]]
  }
  else {
    conservativeIFcalculation[["ic0"]]+getInfluenceFunction.AUC.censoring.term(time = time,
                                                                               event = event,
                                                                               t = t,
                                                                               IFcalculationList = conservativeIFcalculation,
                                                                               MC = MC,
                                                                               cens.model = cens.model, 
                                                                               Wt = Wt,
                                                                               auc = auc, 
                                                                               nth.times = nth.times)
  }
}

getInfluenceCurve.Brier <- function(t,
                                    time,
                                    IC0,
                                    residuals,
                                    IC.G,
                                    cens.model,
                                    nth.times=NULL,
                                    conservative,
                                    event){
  if (is.unsorted(time)){
    stop("Internal error. Time is not sorted in ascending order. ")
  }

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
  else if (cens.model[[1]]=="cox") {## Cox regression
    if (!IC.G$censoring.save.memory){
      ic.weights <- IC.G[[2]][[nth.times]]
      IF.Brier <- IC0+as.numeric(rowSumsCrossprod(as.matrix(residuals),ic.weights,0)) / N
    }
    else {
      wdata <- IC.G[[3]]
      fit <- IC.G[[2]]
      TiMinus <- IC.G[[4]]
      res1 <- (time <= t & event != 0)*residuals
      res2 <- (time > t)*residuals
      IF.Brier <- IC0 + (predictCoxWeights(fit,newdata = wdata,times = TiMinus,diag=TRUE,weights=res1,isBeforeTau = TRUE, tau=t) + predictCoxWeights(fit,newdata = wdata,times = t,diag=FALSE,weights=res2,isBeforeTau = FALSE, tau=t))/N
    }
    IF.Brier
  }else if (cens.model[[1]] == "KaplanMeier"){
    IC0 + getInfluenceFunctionBrierKMCensoringTerm(t,time,residuals,event)
  }
  else {
    stop("Non conservative options with cens.model not being a Cox model or KaplanMeier are not implemented. ")
  }
}
