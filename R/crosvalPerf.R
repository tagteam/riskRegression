# Function to calculate cross-validation performance

crossvalPerf.loob.AUC <- function(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,conservative,cens.model) {
  # initializing output
  if (response.type=="binary")
    auc.loob <- data.table(model=mlevs)
  else
    auc.loob <- data.table(expand.grid(times=times,model=mlevs))
  auc.loob[,AUC:=as.numeric(NA)]
  ## for each pair of individuals sum the concordance of the bootstraps where *both* individuals are out-of-bag
  ## divide by number of times the pair is out-of-bag later
  if (se.fit==TRUE){
    aucDT <- NULL
  }
  ## preparation of outcome status at time horizon(s)
  if (response.type=="binary"){
    NT <- 1
  }
  for (s in 1:NT){
    t <- times[s]
    if (response.type=="binary"){
      ## the following indices have to be logical!!!
      cases.index <- Response[,ReSpOnSe==1]
      controls.index <- !cases.index
      cc.status <- factor(cases.index,levels=c(TRUE,FALSE),labels=c("case","control"))
    }else{
      if (response.type=="survival"){
        ## event of interest before times
        ## the following indices have to be logical!!!
        cases.index <- Response[,time<=t & status==1]
        controls.index <- Response[,time>t]
        cc.status <- factor(rep("censored",N),levels=c("censored","case","control"))
        cc.status[cases.index] <- "case"
        cc.status[controls.index] <- "control"
      }
      else{ ## competing risks
        ## event of interest before times
        ## the following indices have to be logical!!!
        cases.index <- Response[,time<=t & status==1]
        controls.index <- Response[,time>t | status==2]
        cc.status <- factor(rep("censored",N),levels=c("censored","case","control"))
        cc.status[cases.index] <- "case"
        cc.status[controls.index] <- "control"
      }
    }
    # censoring weights
    if (cens.type=="rightCensored"){
      ## IPCW
      weights.cases <- cases.index/Weights$IPCW.subject.times
      if (Weights$method=="marginal"){
        weights.controls <- controls.index/Weights$IPCW.times[s]
      }else{
        weights.controls <- controls.index/Weights$IPCW.times[,s]
      }
      weightMatrix <- outer(weights.cases[cases.index], weights.controls[controls.index], "*")
    }else{ ## uncensored
      weights.cases <- cases.index/1
      weights.controls <- controls.index/1
      weightMatrix <- outer(weights.cases[cases.index], weights.controls[controls.index], "*")
    }
    Phi <- (1/N^2)*sum(weights.cases[cases.index])*sum(weights.controls[controls.index])
    which.cases <- (1:N)[cases.index]
    which.controls <- (1:N)[controls.index]
    for (mod in mlevs){
      Ib <- matrix(0, sum(cases.index), sum(controls.index))
      auc <- matrix(0, sum(cases.index), sum(controls.index))
      if (split.method$internal.name=="crossval"){
        warning("Cannot yet calculate AUC in this case. Use split.method 'loob' or 'bootcv' instead.")
        for (u in 1:B){## cannot use b as running index because b==b does not work in data.table
          riskset <- data.table::data.table(ID=1:N,casecontrol=cc.status,oob=oob)
          data.table::setkey(riskset,ID)
          if (response.type=="binary"){
            oob.risk <- DT.B[model==mod&b==u,data.table::data.table(ID,risk)]
          }else{
            oob.risk <- DT.B[model==mod&times==t&b==u,data.table::data.table(ID,risk)]
          }
          data.table::setkey(oob.risk,ID)
          riskset <- oob.risk[riskset]
          riskset[is.na(risk),risk:=-9]
          Ib.ij <- outer((cases.index*oob)[which.cases],(controls.index*oob)[which.controls],"*")
          auc.ij <- AUCijFun(riskset[casecontrol=="case",risk],riskset[casecontrol=="control",risk])*Ib.ij
          auc <- auc+auc.ij
          auc <- (auc*weightMatrix)/Ib
        }
      }else{
        for (u in 1:B){## cannot use b as running index because b==b does not work in data.table
          ## test <- DT.B[model==mod&times==t&b==u]
          # when B is too low it may happen that some subjects are never oob
          oob <- match(1:N,unique(split.method$index[,u]),nomatch=0)==0
          ## to use the cpp function AUCijFun we
          ## need a vector of length equal to the number of cases (which.cases) for the current time point
          ## which has arbitrary values in places where subjects are inbag and the predicted risks
          ## for those out-of-bag. need another vector for controls.
          riskset <- data.table::data.table(ID=1:N,casecontrol=cc.status,oob=oob)
          data.table::setkey(riskset,ID)
          if (response.type=="binary"){
            oob.risk <- DT.B[model==mod&b==u,data.table::data.table(ID,risk)]
          }else{
            oob.risk <- DT.B[model==mod&times==t&b==u,data.table::data.table(ID,risk)]
          }
          data.table::setkey(oob.risk,ID)
          riskset <- oob.risk[riskset]
          riskset[is.na(risk),risk:=-9]
          Ib.ij <- outer((cases.index*oob)[which.cases],(controls.index*oob)[which.controls],"*")
          auc.ij <- AUCijFun(riskset[casecontrol=="case",risk],riskset[casecontrol=="control",risk])*Ib.ij
          ## Ib.ij is 1 when the pair out of bag
          ## print(head(oob))
          ## print(auc.ij[1:5,1:5])
          auc <- auc+auc.ij
          Ib <- Ib + Ib.ij
        }
        auc <- (auc*weightMatrix)/Ib
      }
      # FIXME: why are there NA's?
      auc[is.na(auc)] <- 0
      ## Leave-one-pair-out bootstrap estimate of AUC
      aucLPO <- (1/N^2)*sum(colSums(auc))*(1/Phi)
      if (is.null(t)){
        auc.loob[model==mod,AUC:=aucLPO]
      }else{
        auc.loob[times==t&model==mod,AUC:=aucLPO]
      }
      if (se.fit==1L){
        ## ## First part of influence function
        ic0Case <- rowSums(auc)
        ic0Control <- colSums(auc)
        ic0 <- (1/(Phi*N))*c(ic0Case, ic0Control)-2*aucLPO
        id.cases <- data[["ID"]][cc.status=="case"]
        id.controls <- data[["ID"]][cc.status=="control"]
        id.censored <- data[["ID"]][cc.status=="censored"]
        if (is.null(t)){
          this.aucDT <- data.table(model=mod,ID = c(id.cases,id.controls), IF.AUC0 = ic0)
        }else{
          this.aucDT <- data.table(model=mod,times=t,ID = c(id.cases,id.controls,id.censored), IF.AUC0 = c(ic0, rep(-2*aucLPO,length(id.censored))))
        }
        aucDT <- rbindlist(list(aucDT,this.aucDT),use.names=TRUE,fill=TRUE)
        if (response.type[[1]]=="binary" || cens.type[[1]]=="uncensored"){
          icPhi <- (aucLPO/Phi)*((weights.cases-(1/N)*sum(weights.cases))*(1/N)*sum(weights.controls)+(weights.controls-(1/N)*sum(weights.controls))*(1/N)*sum(weights.controls))-2*aucLPO
          if (is.null(t)){
            data.table::setkey(aucDT,model,ID)
            aucDT[model==mod,IF.AUC:=IF.AUC0-icPhi]
            auc.loob[model==mod,se:=sd(aucDT[["IF.AUC"]])/sqrt(N)]
          }else{
            data.table::setkey(aucDT,model,times,ID)
            aucDT[times==t&model==mod,IF.AUC:=IF.AUC0-icPhi]
            auc.loob[times==t&model==mod,se:=sd(aucDT[["IF.AUC"]])/sqrt(N)]
          }
        }else{
          ic.weights <- matrix(0,N,N)
          if (cens.type[1]=="rightCensored" && (conservative[1]==FALSE)) {
            ## ## Influence function for G - i.e. censoring survival distribution
            if  (cens.model=="cox"){
              k=0 ## counts subject-times with event before t
              for (i in 1:N){
                if (i %in% id.cases){
                  k=k+1
                  ic.weights[i,] <- Weights$IC$IC.subject[i,k,]/(Weights$IPCW.subject.times[i])
                }else{
                  if (i %in% id.controls){ ## min(T,C)>t
                    ic.weights[i,] <- Weights$IC$IC.times[i,s,]/(Weights$IPCW.times[i,s])
                  }
                }
              }
            }else{
              k=0 ## counts subject-times with event before t
              for (i in 1:N){
                if (i %in% id.cases){
                  ## FIXME: need to check IC
                  pos.i <- sindex(jump.times=unique(data[["time"]]),eval.times=data[["time"]][i])
                  ic.weights[i,] <- Weights$IC[pos.i,]/(Weights$IPCW.subject.times[i])
                }else{
                  if (i %in% id.controls){ ## min(T,C)>t
                    ic.weights[i,] <- Weights$IC[s,]/(Weights$IPCW.times[s])
                  }
                }
              }
            }
            ## ## Part of influence function related to Weights
            ic.weightsCase <- as.numeric(rowSumsCrossprod(as.matrix(rowSums(auc)), ic.weights[which.cases,], 0))
            ic.weightsControl <- as.numeric(rowSumsCrossprod(as.matrix(colSums(auc)), ic.weights[which.controls,], 0))
            ic.weightsCC <- (1/(Phi*N^2))*(ic.weightsCase+ic.weightsControl)
          }
          ## ## Part of influence function related to Phi
          ## icPhiCase <- colMeans(ic.weights[which.cases,])
          icPhiCase <- as.numeric(rowSumsCrossprod(as.matrix(weights.cases[which.cases]),ic.weights[which.cases,],0))
          icPhiControl <- as.numeric(rowSumsCrossprod(as.matrix(weights.controls[which.controls]),ic.weights[which.controls,],0))
          icPhi <- (aucLPO/Phi)*((weights.cases-(1/N)*icPhiCase)*(1/N)*sum(weights.controls)+(weights.controls-(1/N)*icPhiControl)*(1/N)*sum(weights.cases)) - 2*aucLPO
          ## ## Combine all parts of influence function
          ## ic1 <- data.table(ID=data[["ID"]], "ic.weightsCC" = ic.weightsCC, "icPhi" = icPhi)
          data.table::setkey(aucDT,model,times,ID)
          if(conservative==TRUE){
            aucDT[model==mod&times==t, IF.AUC:=IF.AUC0-icPhi]
          }else{
            aucDT[model==mod&times==t, IF.AUC:=IF.AUC0-ic.weightsCC-icPhi]
            aucDT[model==mod&times==t, IF.AUC.conservative:=IF.AUC0-icPhi]
          }
          aucDT[,IF.AUC0:=NULL]
          auc.loob[model==mod&times==t,se:= sd(aucDT[model==mod&times==t,IF.AUC])/sqrt(N)]
          auc.loob[model==mod&times==t,se.conservative:=sd(aucDT[model==mod&times==t,IF.AUC.conservative])/sqrt(N)]
          ## testSE <- sqweightsrt(sum(ic[["ic"]]^2))/N
        }
      }
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

crossvalPerf.loob.Brier <- function(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID){
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
    if (cens.type[1]=="rightCensored" && conservative[1]!=TRUE){
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
        score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,
                                       se=sd(residuals)/sqrt(N),
                                       se.conservative=sd(residuals)/sqrt(N)),by=byvars]
      }else{
        #for small values of B, there is the problem that
        #some individuals might be zero times out of the bag
        #this means that DT.B, which should have a number of rows
        # that is a multiple of the
        #amount of observations in the data, does not fulfill this.
        #the calculations in getInfluenceCurve.Brier cannot accomodate this (for now).


        DT.B[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                time=time,
                                                IC0,
                                                residuals=residuals,
                                                WTi=WTi,
                                                Wt=Wt,
                                                IC.G=Weights$IC,
                                                cens.model=cens.model,
                                                nth.times=nth.times[1]),by=byvars]
        score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,
                                       se=sd(IF.Brier)/sqrt(N),
                                       se.conservative=sd(IC0)/sqrt(N)),by=byvars]


      }
    }else{
      ## either conservative == TRUE or binary or uncensored
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

crossvalPerf.loob <- function(m,times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID) {
  if (m=="AUC"){
    return(crossvalPerf.loob.AUC(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,conservative,cens.model))
  }
  if (m=="Brier"){
    return(crossvalPerf.loob.Brier(times,mlevs,se.fit,response.type,NT,Response,cens.type,Weights,split.method,N,B,DT.B,data,dolist,alpha,byvars,mlabels,ipa,keep.residuals,conservative,cens.model,response.dim,ID))
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
