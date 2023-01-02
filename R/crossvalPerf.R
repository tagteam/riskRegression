# Function to calculate cross-validation performance

crossvalPerf.loob.AUC <- function(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,conservative,cens.model) {
  # initializing output
  bfold <- fold <- AUC <- ReSpOnSe <- status <- ID <- model <- b <- risk <- casecontrol <- IF.AUC <- IF.AUC0 <- se <- IF.AUC.conservative <- se.conservative <- lower <- upper <- NF <- reference  <-  event <- status0 <- NULL
  if (response.type=="binary") {
    auc.loob <- data.table(expand.grid(times=0,model=mlevs)) #add times to auc.loob; now we can write less code for the same thing!
    times <- 0
  }
  auc.loob <- data.table::data.table(expand.grid(times=times,model=mlevs))
  auc.loob[,AUC:=as.numeric(NA)]
  ## for each pair of individuals sum the concordance of the bootstraps where *both* individuals are out-of-bag
  ## divide by number of times the pair is out-of-bag later
  aucDT <- NULL

  redoSplitIndex <- function(x, N){
    res <- rep(TRUE, N)
    res[unique(x)] <- FALSE
    res
  }
  if (split.method$internal.name == "LeaveOneOutBoot"){
    split.index <- apply(split.method$index,2,function(x) redoSplitIndex(x,N))
  }
  else {
    setorder(DT.B, b, ID)
    DT.B$fold <- rep(c(split.method$index), each=NT*length(mlevs))
    DT.B[, bfold:=as.numeric(interaction(b,fold))] # construct new function
    setorder(DT.B, bfold, ID)
    k <- split.method$k
    split.index <- do.call("cbind",lapply(1:k, function(i) split.method$index==i))
  }
  for (s in 1:NT){
    if (response.type=="binary"){
      t <- 0
      ## the following indices have to be logical!!!
      cases.index <- Response[,ReSpOnSe==1]
      controls.index <- !cases.index
    }else{
      t <- times[s]
      if (response.type == "survival"){
        Response$status0 <- Response$status
      }
      else {
        Response$status0 <- Response$status * Response$event
      }
      cases.index <- Response[,time<=t & status0==1]
      controls.index1 <- Response[,time>t]
      controls.index2 <-  Response[,status0==2 & time <= t]
      controls.index <- controls.index1 | controls.index2
    }
    # censoring weights
    if (cens.type=="rightCensored"){ 
      if (Weights$method=="marginal"){
        Wt <- Weights$IPCW.times[s]
      }else{
        Wt <- Weights$IPCW.times[,s]
      }
      weights <- as.numeric(1/Wt * controls.index1  + 1/Weights$IPCW.subject.times * (cases.index | controls.index2))
    }else{ ## uncensored
      weights <- rep(1,N)
    }
    w.cases <- weights[cases.index]
    w.controls <- weights[controls.index]
    n.cases <- sum(cases.index)
    n.controls <- sum(controls.index)
    ID.case <- data[cases.index,]$ID
    ID.controls <- data[controls.index,]$ID
    muCase <- sum(w.cases)
    muControls <- sum(w.controls)
    Phi <- (1/N^2)*muCase*muControls
    for (mod in mlevs){
      if (mod == 0){
        aucLPO <- 0.5
        auc.loob[times==t&model==mod,AUC:=aucLPO]
        DT.B <- DT.B[model != 0]
      }
      else {
        if (response.type == "binary"){ # no reason for the below trick with k fold, as level one data set is approximately N*B
          DTrisk <- DT.B[model == mod][["risk"]]
        }
        else {
          DTrisk <- DT.B[model == mod & times == t][["risk"]]
        }
        if (split.method$internal.name == "crossval"){
          n.col <- k*B
        } else {
          n.col <- B
        }
        risk.mat <- matrix(NA,nrow=N,ncol=n.col)
        risk.mat[split.index] <- DTrisk
        res <- aucLoobFun(ID.case,ID.controls,risk.mat,split.index,weights)
        ic0Case <- res[["ic0Case"]]
        ic0Control <- res[["ic0Control"]]
        warn <- res[["warn"]]
        nu1tauPm <- (1/N^2)*sum(ic0Case)
        ## Leave-one-pair-out bootstrap estimate of AUC
        aucLPO <- nu1tauPm*(1/Phi)
        # following section should be cleaned up(!!!)
        auc.loob[times==t&model==mod,AUC:=aucLPO]
      }
      if (se.fit==1L){ ## should clean this up when cv has been fixed !
        if (mod == 0){
          aucDT <- data.table::data.table(model=mod,times=t,IF.AUC=rep(0,N))
        }
        else {
          IF.AUC0 <- rep(0,N)
          IF.AUC0[cases.index] <- 1/(Phi*N)*ic0Case
          IF.AUC0[controls.index] <- 1/(Phi*N)*ic0Control
          if (!response.type == "binary" && !cens.type == "uncensored" && !conservative[[1]]){
            IFcalculationList <- list(ic0Case = ic0Case,
                                      ic0Control = ic0Control,
                                      weights = weights, 
                                      muCase = muCase,
                                      muControls = muControls, 
                                      nu = nu1tauPm,
                                      firsthit = sindex(data[["time"]],t)-1,
                                      cases = cases.index,
                                      controls =controls.index,
                                      controls1 = controls.index1,
                                      controls2 = controls.index2)
            icPart <- getInfluenceFunction.AUC.censoring.term(time = data[["time"]],
                                                              event = Response[["status0"]],
                                                              t= t,
                                                              IFcalculationList = IFcalculationList,
                                                              MC = Weights$IC, 
                                                              cens.model = cens.model, 
                                                              Wt = Weights$IPCW.times[s], 
                                                              auc = aucLPO, 
                                                              nth.times = s)
          }
          else {
            icPart <- 0
          }
          icPhi1 <- weights*(aucLPO/Phi)*(cases.index*(1/N)*muControls+controls.index*(1/N)*muCase)
          this.aucDT <- data.table::data.table(model=mod,times=t,IF.AUC=IF.AUC0+icPart-icPhi1)
          aucDT <- rbindlist(list(aucDT,this.aucDT),use.names=TRUE,fill=TRUE)
        }
        auc.loob[model==mod&times==t,se:= sd(aucDT[model==mod&times==t,IF.AUC])/sqrt(N)]
      }
    }
  }
  if (warn){
    warning("Some pairs of subjects are never out of bag at the same time. \n You should increase the number of bootstrap replications (argument 'B') ")
  }
  if (response.type == "binary"){
    # remove times again
    auc.loob[,times:=NULL] 
    if (se.fit){
      aucDT[,times:=NULL]
    }
  }
  if (se.fit==1L){
    auc.loob[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
    auc.loob[,upper:=pmin(1,AUC + qnorm(1-alpha/2)*se)]
  }
  output <- list(score=auc.loob)
  if (length(dolist)>0L){
    if (se.fit==FALSE){
      if (match("times",byvars,nomatch=0)){
        contrasts.AUC <- auc.loob[,getComparisons(data.table(x=AUC,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=FALSE),by=list(times)]
      } else{
        contrasts.AUC <- auc.loob[,getComparisons(data.table(x=AUC,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=FALSE)]
      }
    }else{
      aucDT <- merge(aucDT,auc.loob,by=byvars)
      if (match("times",byvars,nomatch=0)){
        contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                               NF=NF,
                                               N=N,
                                               alpha=alpha,
                                               dolist=dolist,
                                               se.fit=TRUE),by=list(times)]
      } else{
        contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                               NF=NF,
                                               N=N,
                                               alpha=alpha,
                                               dolist=dolist,
                                               se.fit=TRUE)]
      }
    }
    setnames(contrasts.AUC,"delta","delta.AUC")
    output <- c(output,list(contrasts=contrasts.AUC))
  }
  if (!is.null(output$score)){
    output$score[,model:=factor(model,levels=mlevs,mlabels)]
  }
  ## set model and reference in model comparison results
  if (!is.null(output$contrasts)>0){
    output$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
    output$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
  }
  return(output)
}

crossvalPerf.loob.Brier <- function(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID,cause){
  status <- residuals <- risk <- WTi <- event <- ReSpOnSe <- Brier <- IC0 <- nth.times <- IF.Brier <- lower <- se <- upper <- model <- NF <- IPCW <- response <- reference <- status0 <- NULL
  ## sum across bootstrap samples where subject i is out of bag
  if (cens.type=="rightCensored"){
    if (response.type=="survival"){
      ## event of interest before times
      DT.B[time<=times & status==1,residuals:=(1-risk)^2/WTi]
    }
    else{ ## competing risks
      ## event of interest before times
      DT.B[time<=times & status==1 & event==cause,residuals:=(1-risk)^2/WTi]
      ## competing event before times
      DT.B[time<=times & status==1 &event!=cause,residuals:=(0-risk)^2/WTi]
    }
    ## right censored before times
    DT.B[time<=times & status==0,residuals:=0]
    ## no event at times
    DT.B[time>times,residuals:=(risk)^2/Wt]
  }else{
    DT.B[,residuals:=(ReSpOnSe-risk)^2]
  }
  ## for each individual sum the residuals of the bootstraps where this individual is out-of-bag
  ## divide by number of times out-off-bag later
  ## DT.B <- DT.B[,data.table::data.table(residuals=sum(residuals)),by=c(byvars,"ID")]
  DT.B <- DT.B[,data.table::data.table(risk=mean(risk),residuals=sum(residuals)),by=c(byvars,"ID")]
  ## get denominator
  if (split.method$name=="LeaveOneOutBoot"){
    Ib <- split.method$B-tabulate(unlist(apply(split.method$index,2,unique)))
    ## REMOVE ME
    ## Ib <- Ib[order(data$ID)]
    if (any(Ib==0)) {
      warning("Some subjects are never out of bag.\n You should increase the number of bootstrap replications (argument 'B').")
      Ib.include <- Ib!=0
      Ib <- Ib[Ib.include]
      ## don't subset residuals, they are only
      ## available for those with Ib.include==1
      ## DT.B <- DT.B[Ib.include]
    } else Ib.include <- NULL
  }else{## cv-k crossvalidation
    Ib <- rep(split.method$B,N)
  }
  ## Ib is the count of how many times subject i is out of bag
  ## the order of Ib matches the order of ID
  ## within groups defined by times and model (byvars)
  data.table::setkeyv(DT.B,c(byvars,"ID"))
  DT.B[,residuals:=residuals/Ib]
  ## leave-one-out bootstrap estimate
  DT.B[,Brier:=mean(residuals),by=byvars]
  ## standard error via influence function
  if (se.fit==1L){
    ## influence function when censoring model is known or data are uncensored
    DT.B[,IC0:=residuals-Brier]
    ## se.brier <- DT.B[,list(se=sd(IC0, na.rm=TRUE)/sqrt(N)),by=byvars]
    ## DT.B[,Brier:=NULL]
    if (cens.type[1]=="rightCensored" && !conservative){
      ## this is a new DT.B
      DT.B <- cbind(data[,1:response.dim,with=FALSE],DT.B)
      DT.B[,nth.times:=as.numeric(factor(times))]
      WW <- data.table(ID=1:N,WTi=Weights$IPCW.subject.times,key="ID")
      DT.B <- merge(WW,DT.B,by=ID)
      ## DT.B[,WTi:=rep(Weights$IPCW.subject.times,NF+length(nullobject))]
      if (Weights$method=="marginal"){
        Wt <- data.table(times=times,Wt=Weights$IPCW.times)
        ## OBS: many digits in times may cause merge problems
        DT.B <- merge(DT.B,Wt,by=c("times"))
      }else{
        Wt <- data.table(times=rep(times,rep(N,NT)),
                         Wt=c(Weights$IPCW.times),
                         ID=data$ID)
        DT.B <- merge(DT.B,Wt,by=c("ID","times"))
      }
      if (cens.type=="uncensored"){
        DT.B[,IF.Brier:= residuals]
        score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=byvars]
      }else{
        #for small values of B, there is the problem that
        #some individuals might be zero times out of the bag
        #this means that DT.B, which should have a number of rows
        # that is a multiple of the
        #amount of observations in the data, does not fulfill this.
        #the calculations in getInfluenceCurve.Brier cannot accomodate this (for now).
        if (response.type == "survival"){
          DT.B[,status0:=status]
        }
        else {
          DT.B[,status0:=status*event]
        }
        DT.B[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                time=time,
                                                IC0,
                                                residuals=residuals,
                                                WTi=WTi,
                                                Wt=Wt,
                                                IC.G=Weights$IC,
                                                cens.model=cens.model,
                                                conservative = conservative,
                                                nth.times=nth.times[1],
                                                event = status0),by=list(model,times)]
        score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,
                                       se=sd(IF.Brier)/sqrt(N)),by=byvars]
      }
    }else{
      ## either binary or uncensored
      score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,se=sd(IC0)/sqrt(N)),
                         by=byvars]
      setnames(DT.B,"IC0","IF.Brier")
    }
    score.loob[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
    score.loob[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
  } else{
    score.loob <- DT.B[,list(Brier=sum(residuals)/N), by=byvars]
  }
  if (ipa==TRUE){
    if (response.type=="binary")
      score.loob[,IPA:=1-Brier/Brier[model==0]]
    else
      score.loob[,IPA:=1-Brier/Brier[model==0],by=times]
  }
  data.table::setkeyv(score.loob,byvars)
  ## data.table::setkey(DT.B,model,times)
  ## DT.B <- DT.B[score.loob]
  if (length(dolist)>0L){
    if (se.fit==FALSE){
      if (match("times",byvars,nomatch=0))
        contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model),
                                                NF=NF,
                                                N=N,
                                                alpha=alpha,
                                                dolist=dolist,
                                                se.fit=FALSE),by=list(times)]
      else
        contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model),
                                                NF=NF,
                                                N=N,
                                                alpha=alpha,
                                                dolist=dolist,
                                                se.fit=FALSE)]
    } else{
      if (match("times",byvars,nomatch=0))
        contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                NF=NF,
                                                N=N,
                                                alpha=alpha,
                                                dolist=dolist,
                                                se.fit=TRUE),by=list(times)]
      else
        contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                NF=NF,
                                                N=N,
                                                alpha=alpha,
                                                dolist=dolist,
                                                se.fit=TRUE)]
    }
    setnames(contrasts.Brier,"delta","delta.Brier")
    output <- list(score=score.loob,contrasts=contrasts.Brier)
  } else{
    output <- list(score=score.loob)
  }
  if (keep.residuals) {
    DT.B[,model:=factor(model,levels=mlevs,mlabels)]
    if (all(c("Wt","WTi")%in%names(DT.B))){
      DT.B[,IPCW:=1/WTi]
      DT.B[time>=times,IPCW:=1/Wt]
      DT.B[time<times & status==0,IPCW:=0]
      output <- c(output,list(residuals=DT.B[,c("ID",names(response),"model","times","risk","residuals","IPCW"),with=FALSE]))
    }else{
      output <- c(output,list(residuals=DT.B[,c("ID",names(response),"model","times","risk","residuals"),with=FALSE]))
    }
  }
  if (!is.null(output$score)){
    output$score[,model:=factor(model,levels=mlevs,mlabels)]
    if (response.type%in%c("survival","competing.risks"))
      setkey(output$score,model,times)
    else
      setkey(output$score,model)
  }
  ## set model and reference in model comparison results
  if (!is.null(output$contrasts)>0){
    output$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
    output$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
    if (response.type%in%c("survival","competing.risks"))
      setkey(output$score,model,times)
    else
      setkey(output$score,model)
  }
  return(output)
}

crossvalPerf.loob <- function(m,times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID,cause) {
  if (m=="AUC"){
    return(crossvalPerf.loob.AUC(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,conservative,cens.model))
  }
  if (m=="Brier"){
    return(crossvalPerf.loob.Brier(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID,cause))
  }
}

crossvalPerf.bootcv <- function(m,crossval,se.fit,multi.split.test,keep.cv,byvars,alpha){
  ## score
  if (length(crossval[[1]][[m]]$score)>0){
    cv.score <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$score}))
    if (se.fit==TRUE){
      bootcv.score <- cv.score[,data.table::data.table(mean(.SD[[m]],na.rm=TRUE),
                                                       se=sd(.SD[[m]],na.rm=TRUE),
                                                       lower=quantile(.SD[[m]],alpha/2,na.rm=TRUE),
                                                       upper=quantile(.SD[[m]],(1-alpha/2),na.rm=TRUE)),by=byvars,.SDcols=m]
      data.table::setnames(bootcv.score,c(byvars,m,"se","lower","upper"))
    }else{
      bootcv.score <- cv.score[,data.table::data.table(mean(.SD[[m]],na.rm=TRUE)),by=byvars,.SDcols=m]
      data.table::setnames(bootcv.score,c(byvars,m))
    }
  }else{
    cv.score <- NULL
    bootcv.score <- NULL
  }
  ## contrasts and multi-split test
  if (length(crossval[[1]][[m]]$contrasts)>0){
    cv.contrasts <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$contrasts}))
    delta.m <- paste0("delta.",m)
    bootcv.contrasts <- switch(as.character(se.fit+3*multi.split.test),
                               "4"={
                                 cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                      lower=quantile(.SD[[delta.m]],alpha/2,na.rm=TRUE),
                                                                      upper=quantile(.SD[[delta.m]],(1-alpha/2),na.rm=TRUE),
                                                                      p=median(.SD[["p"]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m,"p")]
                               },
                               "1"={cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                         lower=quantile(.SD[[delta.m]],alpha/2,na.rm=TRUE),
                                                                         upper=quantile(.SD[[delta.m]],(1-alpha/2),na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m)]
                               },
                               "3"={cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                         p=median(.SD[["p"]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m,"p")]
                               },
                               "0"={
                                 bootcv.contrasts <- cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=delta.m]
                               })
    data.table::setnames(bootcv.contrasts,"V1",delta.m)
  }else{
    cv.contrasts <- NULL
    bootcv.contrasts <- NULL
  }
  out <- list(score=bootcv.score,contrasts=bootcv.contrasts)
  if (keep.cv)
    out <- c(out,list(cv.score=cv.score,cv.contrasts=cv.contrasts))
  out
}
